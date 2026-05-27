# ---------------------------------------------------------------------------
# Target entity type
# ---------------------------------------------------------------------------

"""
    BobbingTarget

Parameters of a vertically oscillating ("bobbing") target.

Fields:
- `y_fixed`: fixed lateral position (m).
- `z_center`: mean altitude (m).
- `z_amplitude`: oscillation amplitude (m).
- `omega`: angular frequency (rad/s).

Call `trajectory(bt, t_range, h, nx, nu)` to materialise the stalk sequence
over any integer time interval.
"""
struct BobbingTarget
    y_fixed::Float64
    z_center::Float64
    z_amplitude::Float64
    omega::Float64
end

"""
    trajectory(bt::BobbingTarget, t_range, h, nx, nu, y_idx, z_idx, zdot_idx)

Return a `Vector` of length `length(t_range)` stalk vectors for `bt`.

Each stalk has dimension `nx + nu`.  The lateral state component at `y_idx`
is set to `bt.y_fixed`; the altitude component at `z_idx` follows a sinusoid
centred at `bt.z_center` with amplitude `bt.z_amplitude` and angular frequency
`bt.omega`; the altitude-rate component at `zdot_idx` is set to the analytic
derivative.  All other components are zero.

`t_range` may be any integer range, enabling evaluation over an arbitrary
sub-interval without recomputing from scratch.
"""
function trajectory(
    bt::BobbingTarget,
    t_range::AbstractVector{<:Integer},
    h::Float64,
    nx::Int,
    nu::Int,
    y_idx::Int,
    z_idx::Int,
    zdot_idx::Int,
)
    @argcheck 1 <= y_idx    <= nx "y index out of range 1:$nx"
    @argcheck 1 <= z_idx    <= nx "z index out of range 1:$nx"
    @argcheck 1 <= zdot_idx <= nx "zdot index out of range 1:$nx"
    return map(t_range) do t
        t_s    = t * h
        z_t    = bt.z_center + bt.z_amplitude * sin(bt.omega * t_s)
        zdot_t = bt.z_amplitude * bt.omega    * cos(bt.omega * t_s)
        x = zeros(nx)
        x[y_idx]    = bt.y_fixed
        x[z_idx]    = z_t
        x[zdot_idx] = zdot_t
        vcat(x, zeros(nu))
    end
end

# ---------------------------------------------------------------------------
# Reference trajectory
# ---------------------------------------------------------------------------

"""
    generate_reference_trajectory(x0, xk, k, Ad, Bd, nx, nu)

Compute the minimum-energy trajectory from `x0` to `xk` in `k` steps by
harmonic extension on a single-agent dynamics sheaf.

Each stalk has dimension `nx + nu`.  Interior control components are in
general nonzero—they are the minimum-energy inputs satisfying the linear
dynamics.  This function is designed to produce fully-pinned target boundary
data; the control components of the returned stalks do not affect the
multi-agent solution once the targets are pinned as boundary vertices.
"""
function generate_reference_trajectory(x0, xk, k, Ad, Bd, nx, nu)
    ref_sheaf = EuclideanSheaf{Float64}(fill(nx + nu, k + 1))
    now_map  = hcat(Matrix(Ad), Matrix(Bd))
    next_map = hcat(Matrix{Float64}(I, nx, nx), zeros(nx, nu))
    for t in 1:k
        add_sheaf_edge!(ref_sheaf, t, t + 1, now_map, next_map)
    end
    bnd = Dict{Int,Vector{Float64}}(
        1     => vcat(Vector{Float64}(x0), zeros(nu)),
        k + 1 => vcat(Vector{Float64}(xk), zeros(nu)),
    )
    z_ref, _ = harmonic_extension(ref_sheaf, bnd)
    return [Array(z_ref[Block(t)]) for t in 1:k+1]
end