# ===========================================================================
# Nonlinear quadrotor on SE(3), with the geometric tracking controller of
#
#   T. Lee, M. Leok and N. H. McClamroch, "Geometric tracking control of a
#   quadrotor UAV on SE(3)," CDC 2010 (arXiv:1003.2005).
#
# This is an implementation of that paper's mathematics, informed by the
# reference implementations at fdcl-gwu/uav_geometric_control (MATLAB) and
# fdcl-gwu/uav_simulator (Python), both MIT-licensed. It is not a port of
# either codebase: no structure, no file layout and no line-for-line
# correspondence is intended, and the frame convention deliberately differs
# (see "Frame convention" below).
#
# WHY THIS CONTROLLER. The layered architecture requires of each agent an ISS
# tracking certificate: a V with a quadratic sandwich and an exponential
# dissipation rate. Lee et al. exhibit exactly that, on an explicit sublevel
# set. So the assumption becomes *checkable* for a real underactuated vehicle
# rather than postulated -- see [`se3_certificate`](@ref), which returns the
# (alpha1, alpha2, c3) triple the composite analysis consumes.
#
# The hover-linearized [`QuadrotorDynamics`](@ref) and [`LQRController`](@ref)
# remain fully supported and are untouched. They are the debugging and
# baselining path: the linear model is the one whose closed loop can be
# checked against a Riccati solution, and any disagreement between the two
# models near hover is a bug in this file.
# ===========================================================================

# ---------------------------------------------------------------------------
# SO(3) utilities
# ---------------------------------------------------------------------------

"""
    hat_so3(w) -> Matrix{Float64}

The isomorphism `R^3 -> so(3)` with `hat_so3(a) * b == cross(a, b)`.
"""
function hat_so3(w::AbstractVector)
    @argcheck length(w) == 3
    return [0.0 -w[3] w[2];
            w[3] 0.0 -w[1];
            -w[2] w[1] 0.0]
end

"""
    vee_so3(S) -> Vector{Float64}

Inverse of [`hat_so3`](@ref). The skew part is taken as given; no symmetry
check is performed, because every caller here feeds it a difference of the
form `A - A'` which is skew by construction.
"""
function vee_so3(S::AbstractMatrix)
    @argcheck size(S) == (3, 3)
    return [S[3, 2], S[1, 3], S[2, 1]]
end

"""
    expm_so3(w) -> Matrix{Float64}

Rodrigues' formula for `exp(hat_so3(w))`. The small-angle branch uses the
Taylor expansions of `sin(t)/t` and `(1 - cos t)/t^2`, which are what keep the
result orthogonal as `norm(w) -> 0`; the closed forms lose all their
significant digits there.
"""
function expm_so3(w::AbstractVector)
    @argcheck length(w) == 3
    theta = norm(w)
    W = hat_so3(w)
    if theta < 1e-8
        # sin(t)/t = 1 - t^2/6, (1 - cos t)/t^2 = 1/2 - t^2/24
        a = 1.0 - theta^2 / 6
        b = 0.5 - theta^2 / 24
    else
        a = sin(theta) / theta
        b = (1 - cos(theta)) / theta^2
    end
    return I + a * W + b * (W * W)
end

"""
    project_SO3(R) -> Matrix{Float64}

Nearest rotation matrix to `R` in Frobenius norm, by SVD. Integrating `Rdot =
R * hat(Omega)` with any classical Runge-Kutta scheme leaves `R` off the
manifold at `O(dt^5)` per step, and the drift accumulates; this projection is
applied once per step to hold it there.
"""
function project_SO3(R::AbstractMatrix)
    @argcheck size(R) == (3, 3)
    F = svd(R)
    Rp = F.U * F.Vt
    if det(Rp) < 0
        # Reflect the least-significant singular direction rather than the
        # first, so a matrix that is already nearly a rotation is perturbed as
        # little as possible.
        S = Diagonal([1.0, 1.0, -1.0])
        Rp = F.U * S * F.Vt
    end
    return Rp
end

"""
    deriv_unit_vector(A, Adot, Addot) -> (q, qdot, qddot)

Normalize `A` and differentiate the normalization, given the derivatives of
`A` itself. Used to turn the commanded force vector and its derivatives into
the desired body-`z` axis and *its* derivatives, which is what supplies the
feedforward angular velocity and acceleration to the attitude loop.

Throws if `norm(A)` is below `tol`: a vanishing commanded force means free
fall, the desired attitude is genuinely undefined there, and silently
returning some arbitrary axis would hide a control design error rather than
surface it.
"""
function deriv_unit_vector(A::AbstractVector, Adot::AbstractVector, Addot::AbstractVector;
                           tol::Float64 = 1e-9)
    nA = norm(A)
    nA > tol || throw(ArgumentError(
        "commanded force vector has norm $nA < $tol: desired attitude is undefined " *
        "(the vehicle is being asked to free-fall). Check the position gains and the " *
        "reference acceleration."))

    q = A / nA
    ndot = dot(A, Adot) / nA
    qdot = Adot / nA - A * (ndot / nA^2)

    nddot = (dot(Adot, Adot) + dot(A, Addot)) / nA - ndot^2 / nA
    qddot = Addot / nA - Adot * (2ndot / nA^2) - A * (nddot / nA^2) +
            A * (2 * ndot^2 / nA^3)

    return (q, qdot, qddot)
end

# ---------------------------------------------------------------------------
# Dynamics
# ---------------------------------------------------------------------------

"""
    SE3QuadrotorDynamics <: AbstractAgentDynamics

Full nonlinear rigid-body quadrotor on `SE(3)`, the plant of Lee et al. (2010):

```
xdot = v
m vdot = m * gravity + f * R * e3
Rdot = R * hat(Omega)
J Omegadot + Omega x J Omega = M
```

with inputs `f` (scalar thrust along the body `z` axis) and `M` (body moment).

# Frame convention

**This module is `z`-up.** `gravity` defaults to `[0, 0, -g]` and thrust acts
along `+R e3`, so at hover `R = I` and `f = m g`. Lee's paper and both FDCL
reference implementations use NED, where `e3` points *down*, gravity is `+g e3`
and the thrust term carries a minus sign. The two are related by a global sign
flip of the vertical axis and the control law below is Lee's, transcribed into
this convention; see [`geometric_control`](@ref) for the term-by-term
correspondence. The choice is deliberate: everything else in this package --
formations, harmonic references, hover altitudes -- is `z`-up, and converting
at the stalk boundary instead would put a frame flip in the hot path of every
experiment, which is where sign errors go to hide.

# State layout

An 18-vector, `[x(1:3); v(4:6); vec(R)(7:15); Omega(16:18)]`, with `R` stored
column-major (`reshape(x[7:15], 3, 3)`). A flat vector rather than a struct so
that the state can travel through the same rollout, logging and filtering
machinery as the linearized models.

# Fields
- `g`: gravitational acceleration magnitude, used to build the default `gravity`.
- `m`: mass.
- `J`: inertia matrix in the body frame, symmetric positive definite.
- `gravity`: the gravity vector itself.
"""
Base.@kwdef struct SE3QuadrotorDynamics <: AbstractAgentDynamics
    g::Float64 = 9.81
    m::Float64 = 0.5
    J::Matrix{Float64} = Matrix(Diagonal([0.01, 0.01, 0.02]))
    gravity::Vector{Float64} = [0.0, 0.0, -9.81]
end

state_dim(::SE3QuadrotorDynamics) = 18
position_indices(::SE3QuadrotorDynamics) = 1:3
velocity_indices(::SE3QuadrotorDynamics) = 4:6

"""
    attitude_indices(dyn::SE3QuadrotorDynamics) -> 7:15
    angular_velocity_indices(dyn::SE3QuadrotorDynamics) -> 16:18
"""
attitude_indices(::SE3QuadrotorDynamics) = 7:15
angular_velocity_indices(::SE3QuadrotorDynamics) = 16:18

"""
    se3_unpack(x) -> (p, v, R, Omega)

Split the 18-vector state. `R` is a view-free copy reshaped to 3x3.
"""
function se3_unpack(x::AbstractVector)
    @argcheck length(x) == 18
    return (x[1:3], x[4:6], reshape(collect(x[7:15]), 3, 3), x[16:18])
end

"""
    se3_pack(p, v, R, Omega) -> Vector{Float64}
"""
function se3_pack(p::AbstractVector, v::AbstractVector, R::AbstractMatrix, Omega::AbstractVector)
    @argcheck length(p) == 3
    @argcheck length(v) == 3
    @argcheck size(R) == (3, 3)
    @argcheck length(Omega) == 3
    return vcat(p, v, vec(R), Omega)
end

"""
    se3_derivative(dyn, x, f, M) -> Vector{Float64}

Time derivative of the 18-vector state under thrust `f` and body moment `M`.
"""
function se3_derivative(dyn::SE3QuadrotorDynamics, x::AbstractVector, f::Real, M::AbstractVector)
    p, v, R, Om = se3_unpack(x)
    e3 = [0.0, 0.0, 1.0]
    vdot = dyn.gravity .+ (f / dyn.m) .* (R * e3)
    Rdot = R * hat_so3(Om)
    Omdot = dyn.J \ (M - cross(Om, dyn.J * Om))
    return vcat(v, vdot, vec(Rdot), Omdot)
end

"""
    se3_step(dyn, x, f, M, dt) -> Vector{Float64}

One step of classical RK4 on [`se3_derivative`](@ref), holding `(f, M)` at
zero order across the interval, followed by [`project_SO3`](@ref) on the
attitude block. Zero-order hold is what a digital autopilot actually does, so
the discretization error this introduces is a property of the implementation
being modelled rather than of the integrator.
"""
function se3_step(dyn::SE3QuadrotorDynamics, x::AbstractVector, f::Real, M::AbstractVector, dt::Real)
    dt > 0 || throw(ArgumentError("dt must be positive"))
    k1 = se3_derivative(dyn, x, f, M)
    k2 = se3_derivative(dyn, x .+ (dt / 2) .* k1, f, M)
    k3 = se3_derivative(dyn, x .+ (dt / 2) .* k2, f, M)
    k4 = se3_derivative(dyn, x .+ dt .* k3, f, M)
    xn = x .+ (dt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    xn[7:15] .= vec(project_SO3(reshape(xn[7:15], 3, 3)))
    return xn
end

"""
    initial_state(dyn::SE3QuadrotorDynamics, position, velocity, acceleration)

Seed a state at `position` with `velocity`, banked into the turn implied by
`acceleration`. The attitude is the one the controller would command there:
body `z` along the total commanded force `m*a - m*gravity`, yaw aligned with
`+x`. This matches the trimming that `initial_state` performs for the
linearized models, so a fleet seeded from the same trajectory starts in the
same physical configuration whichever plant it is flown on.
"""
function initial_state(dyn::SE3QuadrotorDynamics, position::AbstractVector,
                       velocity::AbstractVector, acceleration::AbstractVector)
    @argcheck length(position) == 3
    @argcheck length(velocity) == 3
    a = length(acceleration) >= 3 ? collect(acceleration[1:3]) : zeros(3)
    Fdes = dyn.m .* a .- dyn.m .* dyn.gravity
    R = if norm(Fdes) < 1e-9
        Matrix{Float64}(I, 3, 3)
    else
        b3 = Fdes / norm(Fdes)
        b1d = [1.0, 0.0, 0.0]
        c = cross(b3, b1d)
        if norm(c) < 1e-6
            b1d = [0.0, 1.0, 0.0]
            c = cross(b3, b1d)
        end
        b2 = c / norm(c)
        hcat(cross(b2, b3), b2, b3)
    end
    return se3_pack(position, velocity, R, zeros(3))
end

# ---------------------------------------------------------------------------
# Reference
# ---------------------------------------------------------------------------

"""
    SE3Reference(; x, v, a, jerk, snap, b1, b1dot, b1ddot)

A desired trajectory for [`geometric_control`](@ref).

`x`, `v`, `a` are position, velocity and acceleration and are what the control
law proper needs. `jerk` and `snap` feed *only* the feedforward terms: they
supply the derivatives of the commanded force vector, hence the desired
angular velocity `Omega_c` and acceleration `Omega_c_dot`. Leaving them at
zero is supported and is what a coordination layer that delivers only position
and velocity can do; the cost is that the attitude loop loses its feedforward
and must catch the reference by feedback alone. The stability certificate does
not depend on them.

`b1` is the desired heading direction (not an angle), projected onto the plane
normal to the commanded thrust axis; `b1dot`, `b1ddot` are its derivatives.
"""
Base.@kwdef struct SE3Reference
    x::Vector{Float64} = zeros(3)
    v::Vector{Float64} = zeros(3)
    a::Vector{Float64} = zeros(3)
    jerk::Vector{Float64} = zeros(3)
    snap::Vector{Float64} = zeros(3)
    b1::Vector{Float64} = [1.0, 0.0, 0.0]
    b1dot::Vector{Float64} = zeros(3)
    b1ddot::Vector{Float64} = zeros(3)
end

# ---------------------------------------------------------------------------
# Controller
# ---------------------------------------------------------------------------

"""
    GeometricSE3Controller <: AbstractAgentController

Gains for the geometric tracking law of Lee et al. (2010), equations (19)-(20).

`kx`, `kv` are the position and velocity gains; `kR`, `kOmega` the attitude and
angular-velocity gains. `c1`, `c2` are the Lyapunov cross-term coefficients:
they do not appear in the control law at all, only in the certificate
([`se3_lyapunov`](@ref), [`se3_certificate`](@ref)), where they are what makes
`V` decrease rather than merely be nonincreasing. They are carried on the
controller because the admissible range for them is a function of the gains.

The defaults are stiffer in attitude than the gains shipped with the FDCL
reference implementations. That is not a matter of taste: with `kR` of order
1 and a small quadrotor's inertia, [`se3_certificate`](@ref) certifies an
*empty* region -- the attitude loop is not fast enough for the cross-term
`W12` to be dominated -- so the vehicle flies perfectly well and nothing can be
proved about it. Since a checkable certificate is the reason this controller
is here, the defaults are chosen where one exists.
"""
Base.@kwdef struct GeometricSE3Controller <: AbstractAgentController
    kx::Float64 = 12.0
    kv::Float64 = 8.0
    kR::Float64 = 10.0
    kOmega::Float64 = 2.0
    c1::Float64 = 0.02
    c2::Float64 = 0.06
end

"""
    GeometricSE3Controller(dyn::SE3QuadrotorDynamics; kx, kv, omega_n, zeta)

Gains for a specific vehicle, with the attitude loop specified by the closed-loop
second-order response you want rather than by raw gains.

Linearizing `J * edot_Omega = -kR * eR - kOmega * eOmega` about the identity gives a
second-order system with natural frequency `sqrt(kR / J)` and damping
`kOmega / (2 sqrt(kR J))`, so

```
kR = J * omega_n^2,    kOmega = 2 * zeta * J * omega_n
```

with `J = lambda_max(dyn.J)`.

This constructor exists because of heterogeneity. What sets the attitude loop's
bandwidth is `kR / J`, not `kR`, so a fleet flown on one shared `kR` is detuned in
proportion to how far each vehicle's inertia is from the one the gain was picked
for -- and the heavy vehicles fall out of the certifiable set first, silently.
Specifying `omega_n` instead gives every vehicle the same closed-loop attitude
response and keeps the whole fleet certifiable. See [`se3_certificate`](@ref).
"""
function GeometricSE3Controller(dyn::SE3QuadrotorDynamics;
                                kx::Float64 = 12.0, kv::Float64 = 8.0,
                                omega_n::Float64 = 29.0, zeta::Float64 = 2.9,
                                c1::Float64 = 0.02, c2::Float64 = 0.06)
    Jscalar = maximum(eigvals(Symmetric(dyn.J)))
    return GeometricSE3Controller(kx = kx, kv = kv,
                                  kR = Jscalar * omega_n^2,
                                  kOmega = 2 * zeta * Jscalar * omega_n,
                                  c1 = c1, c2 = c2)
end

"""
    SE3Errors

Tracking errors returned alongside the control input, in the notation of the
paper: `ex`, `ev` translational; `eR`, `eOmega` rotational; `Psi` the attitude
error function `0.5 * tr(I - Rc' * R)`, which is the quantity whose sublevel
sets define the region of attraction. `Rc`, `Omegac` are the commanded
attitude and body rate, retained for logging and for the next step's
finite-difference fallbacks.
"""
struct SE3Errors
    ex::Vector{Float64}
    ev::Vector{Float64}
    eR::Vector{Float64}
    eOmega::Vector{Float64}
    Psi::Float64
    Rc::Matrix{Float64}
    Omegac::Vector{Float64}
end

"""
    geometric_control(ctrl, dyn, state, ref) -> (f, M, err::SE3Errors)

Thrust and body moment for the nonlinear quadrotor `dyn` at `state` tracking
`ref`.

# The law

Translational errors and the commanded force,

```
ex = x - ref.x
ev = v - ref.v
Fdes = -kx * ex - kv * ev - m * gravity + m * ref.a
f = dot(Fdes, R * e3)
b3c = Fdes / norm(Fdes)
```

This is Lee's (19) in `z`-up form. His `A = -kx ex - kv ev - m g e3 + m xddot`
with `e3` down and `f = -dot(A, R e3)`; substituting `gravity = -g e3` and
absorbing the sign gives the above, with `Fdes = -A`. At hover `Fdes` is the
weight and `f = m g`, which is the check worth remembering.

The commanded attitude has `b3c` as its third column and its first column as
close to `ref.b1` as orthogonality allows,

```
b2c = normalize(cross(b3c, ref.b1)),  b1c = cross(b2c, b3c),  Rc = [b1c b2c b3c]
```

and `Omegac = vee(Rc' * Rcdot)`, `Omegacdot = vee(Rc' * Rcddot - hat(Omegac)^2)`,
obtained analytically from the derivatives of `Fdes` -- which is where
`ref.jerk` and `ref.snap` enter, and the only place they enter.

Attitude errors and the moment are Lee's (8), (9), (20),

```
eR = 0.5 * vee(Rc' * R - R' * Rc)
eOmega = Omega - R' * Rc * Omegac
M = -kR * eR - kOmega * eOmega + cross(Omega, J * Omega)
    - J * (hat(Omega) * R' * Rc * Omegac - R' * Rc * Omegacdot)
```

This is the *coupled-yaw* law of the 2010 paper, not the decoupled-yaw law of
Lee (ACC 2019) that the FDCL simulator implements. That is deliberate: the
coupled law is the one whose Lyapunov function and region of attraction are
stated in closed form, and a checkable certificate is the reason this
controller is here.
"""
function geometric_control(ctrl::GeometricSE3Controller, dyn::SE3QuadrotorDynamics,
                           state::AbstractVector, ref::SE3Reference)
    m, J = dyn.m, dyn.J
    e3 = [0.0, 0.0, 1.0]
    x, v, R, Om = se3_unpack(state)

    ex = x .- ref.x
    ev = v .- ref.v

    # Commanded force and thrust, Lee (19) in z-up form.
    Fdes = -ctrl.kx .* ex .- ctrl.kv .* ev .- m .* dyn.gravity .+ m .* ref.a
    b3 = R * e3
    f = dot(Fdes, b3)

    # Derivatives of the commanded force. `ea` is the actual acceleration
    # error and `eb` its derivative; both are the realized quantities, not the
    # commanded ones, which is what makes the feedforward exact rather than
    # approximate.
    a_actual = dyn.gravity .+ (f / m) .* b3
    ea = a_actual .- ref.a
    Fdes_dot = -ctrl.kx .* ev .- ctrl.kv .* ea .+ m .* ref.jerk

    b3_dot = R * hat_so3(Om) * e3
    f_dot = dot(Fdes_dot, b3) + dot(Fdes, b3_dot)
    eb = (f_dot / m) .* b3 .+ (f / m) .* b3_dot .- ref.jerk
    Fdes_ddot = -ctrl.kx .* ea .- ctrl.kv .* eb .+ m .* ref.snap

    b3c, b3c_dot, b3c_ddot = deriv_unit_vector(Fdes, Fdes_dot, Fdes_ddot)

    # Second body axis from the heading direction, and its derivatives.
    A2 = cross(b3c, ref.b1)
    A2_dot = cross(b3c_dot, ref.b1) .+ cross(b3c, ref.b1dot)
    A2_ddot = cross(b3c_ddot, ref.b1) .+ 2 .* cross(b3c_dot, ref.b1dot) .+
              cross(b3c, ref.b1ddot)
    b2c, b2c_dot, b2c_ddot = deriv_unit_vector(A2, A2_dot, A2_ddot)

    b1c = cross(b2c, b3c)
    b1c_dot = cross(b2c_dot, b3c) .+ cross(b2c, b3c_dot)
    b1c_ddot = cross(b2c_ddot, b3c) .+ 2 .* cross(b2c_dot, b3c_dot) .+ cross(b2c, b3c_ddot)

    Rc = hcat(b1c, b2c, b3c)
    Rc_dot = hcat(b1c_dot, b2c_dot, b3c_dot)
    Rc_ddot = hcat(b1c_ddot, b2c_ddot, b3c_ddot)

    Omegac = vee_so3(Rc' * Rc_dot)
    Omegac_dot = vee_so3(Rc' * Rc_ddot .- hat_so3(Omegac)^2)

    # Attitude errors, Lee (6), (8), (9).
    eR = 0.5 .* vee_so3(Rc' * R .- R' * Rc)
    RtRc = R' * Rc
    eOmega = Om .- RtRc * Omegac
    Psi = 0.5 * tr(I - Rc' * R)

    # Moment, Lee (20).
    M = -ctrl.kR .* eR .- ctrl.kOmega .* eOmega .+ cross(Om, J * Om) .-
        J * (hat_so3(Om) * RtRc * Omegac .- RtRc * Omegac_dot)

    return (f, M, SE3Errors(ex, ev, eR, eOmega, Psi, Rc, Omegac))
end

# ---------------------------------------------------------------------------
# Lyapunov certificate
# ---------------------------------------------------------------------------

"""
    se3_lyapunov(ctrl, dyn, err) -> Float64

The Lyapunov function of Lee et al., Proposition 3,

```
V = 0.5 * kx * |ex|^2 + 0.5 * m * |ev|^2 + c1 * dot(ex, ev)
  + 0.5 * dot(eOmega, J * eOmega) + kR * Psi + c2 * dot(eR, eOmega)
```

evaluated on the errors returned by [`geometric_control`](@ref). Positive
definite on the region where the gain conditions of
[`se3_certificate`](@ref) hold; note that it uses `Psi` and not `|eR|^2`, the
two agreeing only to second order.
"""
function se3_lyapunov(ctrl::GeometricSE3Controller, dyn::SE3QuadrotorDynamics, err::SE3Errors)
    V1 = 0.5 * ctrl.kx * dot(err.ex, err.ex) + 0.5 * dyn.m * dot(err.ev, err.ev) +
         ctrl.c1 * dot(err.ex, err.ev)
    V2 = 0.5 * dot(err.eOmega, dyn.J * err.eOmega) + ctrl.kR * err.Psi +
         ctrl.c2 * dot(err.eR, err.eOmega)
    return V1 + V2
end

"""
    SE3Certificate

What Proposition 3 of Lee et al. yields, in the form the layered stability
analysis consumes.

- `alpha1`, `alpha2`: the quadratic sandwich constants,
  `alpha1 * |z|^2 <= V <= alpha2 * |z|^2`, with
  `z = (|ex|, |ev|, |eR|, |eOmega|)`.
- `c3`: the exponential dissipation rate, `Vdot <= -c3 * V`.
- `feasible`: whether every hypothesis actually holds at these gains. **If
  this is false the other three numbers are not a certificate of anything.**
- `psi1`: the attitude sublevel set the certificate is stated on. `V` is a
  valid certificate on `{Psi < psi1}` intersected with the position bound.
- `region_ex_max`: the position-error bound the argument was closed under.
- `violations`: the *binding* conditions that failed. Empty iff `feasible`.
- `notes`: advisory failures of Lee's closed-form sufficient bounds on `c1`
  and `c2`. Those bounds are sufficient for the binding conditions but not
  necessary, so a certificate can be perfectly valid with notes attached; they
  are reported because they say *which* gain is near its limit.

The region is the point of this type. No global certificate exists for this
plant -- with attitude on `SO(3)` and thrust of one sign there is no globally
valid quadratic sandwich -- so any composite result built on it is regional,
and this struct is where that region is carried rather than assumed away.
"""
struct SE3Certificate
    alpha1::Float64
    alpha2::Float64
    c3::Float64
    feasible::Bool
    psi1::Float64
    region_ex_max::Float64
    violations::Vector{String}
    notes::Vector{String}
end

"""
    se3_certificate(ctrl, dyn; psi1, ex_max, ev_max, a_max) -> SE3Certificate

Evaluate Proposition 3 of Lee et al. (2010) at these gains and this vehicle.

The proposition bounds `V` between two quadratic forms and bounds `Vdot` by a
third, in the stacked variable `z = (|ex|, |ev|, |eR|, |eOmega|)`:

```
z' * blockdiag(M11, M21) * z <= V <= z' * blockdiag(M12, M22) * z
Vdot <= -z1' * W1 * z1 + z1' * W12 * z2 - z2' * W2 * z2
```

with, writing `alpha = sqrt(psi1 * (2 - psi1))`, `Jm = lambda_min(J)`,
`JM = lambda_max(J)` and `B` a bound on the commanded force,

```
M11 = 0.5 * [kx  -c1;  -c1  m]        M12 = 0.5 * [kx  c1;  c1  m]
M21 = 0.5 * [2kR -c2;  -c2  Jm]       M22 = 0.5 * [2kR/(2-psi1)  c2;  c2  JM]

W1  = [c1*kx*(1-alpha)/m        -c1*kv*(1+alpha)/(2m);
       -c1*kv*(1+alpha)/(2m)     kv*(1-alpha) - c1]
W2  = [c2*kR/JM                 -c2*kOmega/(2*Jm);
       -c2*kOmega/(2*Jm)         kOmega - c2]
W12 = [c1*B/m  0;  B + kx*ex_max  0]
```

Rather than checking Lee's closed-form sufficient inequalities on `c1` and
`c2` and stopping there, this assembles the full quadratic form

```
W = [W1  -W12/2;  -W12'/2  W2]
```

and tests `W > 0` numerically. That is the condition the closed forms are
sufficient for, so this accepts every gain set they accept and some they
reject conservatively. The returned rate is `c3 = lambda_min(W) /
lambda_max(blockdiag(M12, M22))`, the standard conversion from
`Vdot <= -z'Wz` and `V <= z'Mz`.

# Arguments
- `psi1`: attitude sublevel set, in `(0, 1)`. Lee's Proposition 3 requires
  `psi1 < 1`; the almost-global result of Proposition 4 covers `psi1 < 2` but
  gives attractiveness rather than an exponential estimate, so it is not what
  a composite ISS argument can use.
- `ex_max`: position-error bound defining the region.
- `ev_max`, `a_max`: velocity-error and reference-acceleration bounds, used
  only to form `B = kx*ex_max + kv*ev_max + m*g + m*a_max`, the bound on the
  commanded force that drives the cross-term `W12`.

# Example

```julia
dyn  = SE3QuadrotorDynamics(m = 0.6, J = Matrix(Diagonal([0.012, 0.012, 0.024])))
ctrl = GeometricSE3Controller(kx = 12.0, kv = 8.0, kR = 1.8, kOmega = 0.45,
                              c1 = 1.2, c2 = 0.30)
cert = se3_certificate(ctrl, dyn; psi1 = 0.6, ex_max = 1.0, ev_max = 1.0, a_max = 1.0)
cert.feasible && cert.c3   # the Assumption-1 decay rate for this agent
```
"""
function se3_certificate(ctrl::GeometricSE3Controller, dyn::SE3QuadrotorDynamics;
                         psi1::Float64 = 0.9,
                         ex_max::Float64 = 1.0,
                         ev_max::Float64 = 1.0,
                         a_max::Float64 = 0.0)
    violations = String[]
    0 < psi1 < 1 || push!(violations, "psi1 = $psi1 not in (0,1); Prop. 3 requires psi1 < 1")
    ex_max > 0 || push!(violations, "ex_max must be positive")

    m = dyn.m
    Jm = minimum(eigvals(Symmetric(dyn.J)))
    JM = maximum(eigvals(Symmetric(dyn.J)))
    Jm > 0 || push!(violations, "J must be positive definite")

    alpha = sqrt(clamp(psi1 * (2 - psi1), 0.0, 1.0))
    kx, kv, kR, kOm, c1, c2 = ctrl.kx, ctrl.kv, ctrl.kR, ctrl.kOmega, ctrl.c1, ctrl.c2

    B = kx * ex_max + kv * ev_max + m * norm(dyn.gravity) + m * a_max

    M11 = 0.5 .* [kx -c1; -c1 m]
    M12 = 0.5 .* [kx c1; c1 m]
    M21 = 0.5 .* [2kR -c2; -c2 Jm]
    M22 = 0.5 .* [2kR / (2 - psi1) c2; c2 JM]

    W1 = [c1 * kx * (1 - alpha)/m        -c1 * kv * (1 + alpha)/(2m);
          -c1 * kv * (1 + alpha)/(2m)     kv * (1 - alpha) - c1]
    W2 = [c2 * kR / JM            -c2 * kOm / (2Jm);
          -c2 * kOm / (2Jm)        kOm - c2]
    W12 = [c1 * B / m  0.0;
           B + kx * ex_max  0.0]

    # Lee's closed-form sufficient conditions on c1 and c2. These are
    # sufficient for the quadratic forms below to have the right sign, not
    # necessary, and they are markedly conservative -- so they are recorded as
    # advisory notes and the binding test is positive definiteness of the
    # assembled forms.
    c1_bound = min(kv * (1 - alpha),
                   4m * kx * kv * (1 - alpha)^2 /
                       (kv^2 * (1 + alpha)^2 + 4m * kx * (1 - alpha)),
                   sqrt(kx * m))
    c2_bound = min(kOm,
                   4kOm * kR * Jm^2 / (kOm^2 * JM + 4kR * Jm^2),
                   sqrt(kR * Jm))
    notes = String[]
    c1 < c1_bound || push!(notes, "c1 = $c1 exceeds Lee's sufficient bound $(round(c1_bound, digits=4))")
    c2 < c2_bound || push!(notes, "c2 = $c2 exceeds Lee's sufficient bound $(round(c2_bound, digits=4))")

    lower = _blockdiag2(M11, M21)
    upper = _blockdiag2(M12, M22)
    W = [W1 (-0.5 .* W12); (-0.5 .* W12') W2]

    lam_lower = minimum(eigvals(Symmetric(lower)))
    lam_upper = maximum(eigvals(Symmetric(upper)))
    lam_W = minimum(eigvals(Symmetric(W)))

    lam_lower > 0 || push!(violations, "V is not positive definite: lambda_min = $(round(lam_lower, sigdigits=3))")
    lam_W > 0 || push!(violations, "Vdot bound is not negative definite: lambda_min(W) = $(round(lam_W, sigdigits=3))")

    feasible = isempty(violations)
    c3 = feasible ? lam_W / lam_upper : 0.0

    return SE3Certificate(max(lam_lower, 0.0), lam_upper, c3, feasible,
                          psi1, ex_max, violations, notes)
end

"""Assemble `[A 0; 0 B]` for the 2x2 blocks used by the certificate, with the
rows ordered `(ex, ev, eR, eOmega)`."""
function _blockdiag2(A::AbstractMatrix, B::AbstractMatrix)
    Z = zeros(2, 2)
    return [A[1,1] A[1,2] 0.0 0.0;
            A[2,1] A[2,2] 0.0 0.0;
            0.0 0.0 B[1,1] B[1,2];
            0.0 0.0 B[2,1] B[2,2]]
end

"""
    se3_error_vector(err::SE3Errors) -> Vector{Float64}

The stacked magnitudes `z = (|ex|, |ev|, |eR|, |eOmega|)` that the certificate
of [`se3_certificate`](@ref) is written in. Provided so that a rollout can
check `alpha1*|z|^2 <= V <= alpha2*|z|^2` and `Vdot <= -c3*V` numerically
along a real trajectory rather than trusting the algebra.
"""
function se3_error_vector(err::SE3Errors)
    return [norm(err.ex), norm(err.ev), norm(err.eR), norm(err.eOmega)]
end

"""
    se3_tune_lyapunov(ctrl, dyn; psi1, ex_max, ev_max, a_max, n, widen) -> GeometricSE3Controller

Search for `(c1, c2)` maximizing the certified decay rate `c3` at fixed
control gains, by grid search. Returns a controller with the best pair
substituted; the control law is unchanged, since `c1` and `c2` never enter it.

The box searched is `widen` times Lee's closed-form admissible box. Because
those closed forms are sufficient but not necessary and the binding test in
[`se3_certificate`](@ref) is numerical, `widen > 1` legitimately finds
certificates that the closed forms reject; the returned certificate carries a
`note` saying so. Set `widen = 1` to stay strictly inside the published
bounds.

Returns `ctrl` unchanged if no feasible pair exists in the box; call
[`se3_autotune`](@ref) instead of checking that by hand.

Worth running rather than trusting the defaults: the certified rate is
strongly sensitive to the cross-term coefficients, and the struct defaults are
chosen to be safe across vehicles rather than tight for any one of them.
"""
function se3_tune_lyapunov(ctrl::GeometricSE3Controller, dyn::SE3QuadrotorDynamics;
                           psi1::Float64 = 0.9, ex_max::Float64 = 1.0,
                           ev_max::Float64 = 1.0, a_max::Float64 = 0.0, n::Int = 40,
                           widen::Float64 = 4.0)
    Jm = minimum(eigvals(Symmetric(dyn.J)))
    JM = maximum(eigvals(Symmetric(dyn.J)))
    alpha = sqrt(clamp(psi1 * (2 - psi1), 0.0, 1.0))
    m, kx, kv, kR, kOm = dyn.m, ctrl.kx, ctrl.kv, ctrl.kR, ctrl.kOmega

    c1_hi = widen * min(kv * (1 - alpha),
                        4m * kx * kv * (1 - alpha)^2 /
                            (kv^2 * (1 + alpha)^2 + 4m * kx * (1 - alpha)),
                        sqrt(kx * m))
    c2_hi = widen * min(kOm, 4kOm * kR * Jm^2 / (kOm^2 * JM + 4kR * Jm^2), sqrt(kR * Jm))

    best = ctrl
    best_c3 = -Inf
    for c1 in range(c1_hi / n, c1_hi * (1 - 1 / n); length = n),
        c2 in range(c2_hi / n, c2_hi * (1 - 1 / n); length = n)

        trial = GeometricSE3Controller(; kx, kv, kR, kOmega = kOm, c1, c2)
        cert = se3_certificate(trial, dyn; psi1, ex_max, ev_max, a_max)
        if cert.feasible && cert.c3 > best_c3
            best_c3 = cert.c3
            best = trial
        end
    end
    return best
end

"""
    se3_autotune(ctrl, dyn; kwargs...) -> (GeometricSE3Controller, SE3Certificate)

[`se3_tune_lyapunov`](@ref) together with the certificate it achieved. Prefer
this at call sites: `se3_tune_lyapunov` returns its input unchanged when no
feasible `(c1, c2)` exists anywhere in the search box, and an unchanged
controller is indistinguishable from a successfully tuned one until you look
at the certificate. Here the certificate comes back with it, so
`cert.feasible` cannot be skipped by accident.

Keyword arguments are those of [`se3_tune_lyapunov`](@ref).
"""
function se3_autotune(ctrl::GeometricSE3Controller, dyn::SE3QuadrotorDynamics; kwargs...)
    tuned = se3_tune_lyapunov(ctrl, dyn; kwargs...)
    cert_kwargs = filter(kv -> kv.first in (:psi1, :ex_max, :ev_max, :a_max), kwargs)
    return (tuned, se3_certificate(tuned, dyn; cert_kwargs...))
end

"""
    se3_measured_rate(ctrl, dyn, x0, ref, dt, nsteps) -> (rate, Vs)

Realized exponential decay rate of `V` along an actual rollout from `x0`
holding `ref` fixed: the least-squares slope of `log V` against time, negated.

This exists because the certified `c3` of [`se3_certificate`](@ref) is a
sufficient bound and, for this Lyapunov function, a very loose one -- an order
of magnitude or two below the rate the closed loop actually achieves. Both
numbers are worth having and they answer different questions. `c3` is what a
composite stability argument may assume; this is what the vehicle does.

Returns the rate and the sampled `V` trajectory. The fit uses only samples
where `V` is above `floor_frac` of its initial value, so that the numerical
floor at the end of the run does not flatten the slope.
"""
function se3_measured_rate(ctrl::GeometricSE3Controller, dyn::SE3QuadrotorDynamics,
                           x0::AbstractVector, ref::SE3Reference, dt::Real, nsteps::Int;
                           floor_frac::Float64 = 1e-6)
    x = collect(float.(x0))
    Vs = Float64[]
    for _ in 1:nsteps
        f, M, err = geometric_control(ctrl, dyn, x, ref)
        push!(Vs, se3_lyapunov(ctrl, dyn, err))
        x = se3_step(dyn, x, f, M, dt)
    end
    isempty(Vs) && return (0.0, Vs)

    keep = findall(v -> v > floor_frac * Vs[1] && v > 0, Vs)
    length(keep) < 2 && return (0.0, Vs)

    t = (keep .- 1) .* float(dt)
    y = log.(Vs[keep])
    tbar, ybar = sum(t) / length(t), sum(y) / length(y)
    num = sum((t .- tbar) .* (y .- ybar))
    den = sum((t .- tbar) .^ 2)
    return (den > 0 ? -num / den : 0.0, Vs)
end

# ---------------------------------------------------------------------------
# Agent wrapper: reference filter + geometric controller + nonlinear plant
# ---------------------------------------------------------------------------

"""
    SE3AgentState <: AbstractAgentState

One nonlinear agent in the layered stack: a [`JointTikhonovFilter`](@ref) on
the delivered coordination reference, a [`GeometricSE3Controller`](@ref), and
an [`SE3QuadrotorDynamics`](@ref) plant.

The filter runs on the three *position* coordinates only. That is the whole
point of the layer separation: the coordination layer emits a position
reference, the filter smooths it, and what happens between that reference and
the rotors is the agent's own business. Filtering a rotation matrix
elementwise, which is what reusing the flat-state filter of the linearized
agents would do here, is not a meaningful operation.

# Fields
- `x`: 18-vector plant state.
- `filter`: joint position/velocity reference filter.
- `ctrl`, `dyn`: controller gains and vehicle.
- `b1`: desired heading direction, held fixed by default.
- `feedforward`: when true, the filter is driven by `qstar + eps * qstar_dot`
  rather than by `qstar`, which cancels the moving-reference forcing term
  identically and so removes the filter's own `O(eps)` lag. This is the
  analytic feedforward of the reference-filter proposition, and it is a
  property of the *agent* only in the sense of which signal it is handed --
  the extra right-hand side is solved for upstream, in the same two sweeps of
  the same cached factor.
- `last`: errors from the most recent step, for logging and certificate checks.
"""
mutable struct SE3AgentState <: AbstractAgentState
    x::Vector{Float64}
    filter::JointTikhonovFilter{Float64, Vector{Float64}}
    ctrl::GeometricSE3Controller
    dyn::SE3QuadrotorDynamics
    b1::Vector{Float64}
    feedforward::Bool
    last::Union{Nothing, SE3Errors}
end

"""
    SE3AgentState(x0, dyn, ctrl, eps; b1, r0, rdot0)

Build an agent at plant state `x0` with reference-filter time constant `eps`.
The filter is seeded at `r0` / `rdot0`, defaulting to the agent's own position
and velocity -- which sets the planner error `e(0)` to zero, so the composite
transient shows the *vehicle* converging to the formation rather than the
filter converging to the vehicle. Seeding it elsewhere is supported and is how
the cascade's peaking behaviour can be exercised deliberately.
"""
function SE3AgentState(x0::AbstractVector, dyn::SE3QuadrotorDynamics,
                       ctrl::GeometricSE3Controller, eps::Float64;
                       b1::AbstractVector = [1.0, 0.0, 0.0],
                       feedforward::Bool = false,
                       r0::Union{Nothing, AbstractVector} = nothing,
                       rdot0::Union{Nothing, AbstractVector} = nothing)
    @argcheck length(x0) == 18
    p0 = r0 === nothing ? collect(x0[1:3]) : collect(r0)
    v0 = rdot0 === nothing ? collect(x0[4:6]) : collect(rdot0)
    flt = JointTikhonovFilter(p0, v0; epsilon = eps)
    return SE3AgentState(collect(float.(x0)), flt, ctrl, dyn, collect(float.(b1)),
                         feedforward, nothing)
end

"""
    step_agent!(w::SE3AgentState, qstar, qstar_dot, qstar_ddot, dt) -> (x, ref)

Advance one control step: filter the delivered reference, compute `(f, M)`
from the geometric law, and integrate the nonlinear plant under zero-order
hold.

`qstar` is agent `i`'s block of the harmonic extension. Only the first three
coordinates are used, so a stalk carrying SE(3) homogeneous coordinates -- a
position and a trailing constant 1 -- may be passed unmodified. `qstar_dot`
and `qstar_ddot` are the corresponding reference velocity and acceleration;
`qstar_ddot` is passed through to the controller unfiltered, as the
acceleration feedforward, and may be omitted.

Returns the new plant state and the filtered reference position actually
tracked.
"""
function step_agent!(w::SE3AgentState, qstar::AbstractVector, qstar_dot::AbstractVector,
                     qstar_ddot::AbstractVector, dt::Real)
    q = _first3(qstar)
    qd = _first3(qstar_dot)
    qdd = _first3(qstar_ddot)

    # Analytic feedforward: driving the filter with qstar + eps*qstar_dot makes
    # its error dynamics eps*edot = -e, so the moving-reference term cancels
    # identically instead of being tracked with an O(eps) lag. The velocity
    # channel gets the same treatment one derivative up.
    if w.feedforward
        tikhonov_step!(w.filter, q .+ w.filter.epsilon .* qd,
                       qd .+ w.filter.epsilon .* qdd, dt)
    else
        tikhonov_step!(w.filter, q, qd, dt)
    end
    ref = SE3Reference(x = copy(w.filter.x), v = copy(w.filter.v), a = qdd, b1 = w.b1)

    f, M, err = geometric_control(w.ctrl, w.dyn, w.x, ref)
    w.x = se3_step(w.dyn, w.x, f, M, dt)
    w.last = err
    return (copy(w.x), copy(ref.x))
end

function step_agent!(w::SE3AgentState, qstar::AbstractVector, qstar_dot::AbstractVector, dt::Real)
    return step_agent!(w, qstar, qstar_dot, zeros(3), dt)
end

function step_agent!(w::SE3AgentState, qstar::AbstractVector, dt::Real)
    return step_agent!(w, qstar, zeros(3), zeros(3), dt)
end

"""Take the leading three coordinates of a stalk value, padding with zeros if
it is shorter. Stalks in this package are commonly 4-dimensional (SE(3)
homogeneous coordinates, whose fourth entry is the constant 1); the trailing
entry carries no kinematic information and is dropped here rather than at
every call site."""
function _first3(q::AbstractVector)
    n = length(q)
    n >= 3 && return collect(float.(q[1:3]))
    out = zeros(3)
    out[1:n] .= q
    return out
end
