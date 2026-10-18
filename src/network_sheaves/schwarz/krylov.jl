# Krylov acceleration of Schwarz sweeps for SchwarzMethods: GMRES and BiCGStab,
# preconditioned by a sweep. Included in SchwarzMethods.

export SchwarzGMRES, SchwarzSweepPreconditioner, SchwarzBiCGStab

# ===== GMRES acceleration =====

"""
    SchwarzSweepPreconditioner(dd, sweep, coarse=nothing)

One step of a stationary Schwarz iteration as a linear map: `P * r` starts from
the zero cochain, applies one `sweep` of [`schwarz_step!`](@ref) to
``A e = r``, then the `coarse` correction ([`coarse_correct!`](@ref)) if any,
and returns the glued result ``e = M^{-1} r``. The map is linear in ``r`` for
every sweep, transmission condition and coarse space:

- with [`ParallelSweep`](@ref) and Dirichlet transmission ``M^{-1}`` is
  restricted additive Schwarz (RAS);
- with [`ParallelSweep`](@ref) and [`RobinTransmission`](@ref) it is optimized
  RAS (ORAS, St-Cyr–Gander–Thomas 2007);
- with a coarse space it is the corresponding hybrid two-level method,
  ``M^{-1} = M_1^{-1} + \\Phi A_0^{-1}\\Phi^\\mathsf{T}(I - A M_1^{-1})`` for a
  [`TruncatedPushforwardCoarseSpace`](@ref).

It is not symmetric in general, so it preconditions GMRES ([`SchwarzGMRES`](@ref)).
"""
struct SchwarzSweepPreconditioner{T,D<:SchwarzDecomposition{T},S<:SchwarzSweep,C<:Union{Nothing,AbstractCoarseSpace}}
    decomposition::D
    sweep::S
    coarse::C
    owned::Vector{Vector{Int}}      # positions t in subdomain i with owner(dofs[t]) == i
    rhs::Vector{Vector{T}}          # per-subdomain buffers for the direct RAS apply
    solution::Vector{Vector{T}}
end

function SchwarzSweepPreconditioner(dd::SchwarzDecomposition{T}, sweep::SchwarzSweep, coarse=nothing) where {T}
    owned = Vector{Vector{Int}}(undef, length(dd.locals))
    Threads.@threads for i in eachindex(owned)
        dofs = dd.locals[i].dofs
        owned[i] = [t for t in eachindex(dofs) if dd.ownership.owner[dofs[t]] == i]
    end
    rhs = [zeros(T, length(lp.dofs)) for lp in dd.locals]
    return SchwarzSweepPreconditioner{T,typeof(dd),typeof(sweep),typeof(coarse)}(dd, sweep, coarse, owned,
        rhs, map(similar, rhs))
end

Base.size(P::SchwarzSweepPreconditioner) = size(P.decomposition.A)
Base.size(P::SchwarzSweepPreconditioner, d::Integer) = size(P.decomposition.A, d)
Base.eltype(::SchwarzSweepPreconditioner{T}) where {T} = T

# Restricted additive Schwarz from the zero cochain: the ghost values are zero,
# so each local solve sees r on its dofs only, and gluing keeps the owned
# entries. Same map as the general path below, without building the cochain.
function _ras!(y::AbstractVector, P::SchwarzSweepPreconditioner, r::AbstractVector)
    dd = P.decomposition
    Threads.@threads for i in eachindex(dd.locals)
        lp = dd.locals[i]
        b, z = P.rhs[i], P.solution[i]
        @inbounds for t in eachindex(b)
            b[t] = r[lp.dofs[t]]
        end
        _factor_solve!(z, lp.factor, b)
        @inbounds for t in P.owned[i]
            y[lp.dofs[t]] = z[t]
        end
    end
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, P::SchwarzSweepPreconditioner{T}, r::AbstractVector) where {T}
    dd = P.decomposition
    @argcheck length(y) == length(r) == size(dd.A, 1)
    P.sweep isa ParallelSweep && P.coarse === nothing && return _ras!(y, P, r)
    n = size(dd.A, 1)
    prob = SchwarzProblem(dd, Vector{T}(r), zeros(T, n))
    x = localize(dd, prob.u0)
    schwarz_step!(x, prob, P.sweep)
    P.coarse === nothing || coarse_correct!(x, prob, P.coarse)
    y .= glue(dd, x)
    return y
end

Base.:*(P::SchwarzSweepPreconditioner{T}, r::AbstractVector) where {T} =
    mul!(similar(r, promote_type(T, eltype(r))), P, r)

"""
    SchwarzGMRES(; sweep=ParallelSweep(), coarse=nothing, tol=1e-8, maxiter=1000, memory=100)

GMRES (Krylov.jl, restarted every `memory` iterations) right-preconditioned
with one step of the stationary Schwarz iteration
([`SchwarzSweepPreconditioner`](@ref)). It accelerates any
[`SchwarzIteration`](@ref) with the same sweep and coarse space, including
combinations whose stationary iteration diverges. The standard example is
optimized (Robin) RAS with a coarse level, the two-level ORAS method
(Dolean–Jolivet–Nataf 2015, ch. 5).

Right preconditioning makes GMRES minimize the true residual. Iteration stops
when ``\\|f - A u\\| \\le \\mathrm{tol}\\,\\|f\\|``, the same criterion as
[`SchwarzIteration`](@ref) and [`SchwarzCG`](@ref). The relative residual of
the final iterate is checked and recorded.
"""
Base.@kwdef struct SchwarzGMRES{S<:SchwarzSweep,C<:Union{Nothing,AbstractCoarseSpace}}
    sweep::S = ParallelSweep()
    coarse::C = nothing
    tol::Float64 = 1e-8
    maxiter::Int = 1000
    memory::Int = 100
end

function solve(prob::SchwarzProblem{T}, alg::SchwarzGMRES) where {T}
    @argcheck alg.maxiter >= 1 && alg.memory >= 1
    dd = prob.decomposition
    r0 = prob.rhs - _threaded_mul!(similar(prob.rhs), dd.At, prob.u0)
    initial = _relative_residual(prob, prob.u0)
    if initial <= alg.tol || iszero(norm(r0))
        return SchwarzResult(copy(prob.u0), localize(dd, prob.u0), T[initial], T[0], 0, true)
    end
    # GMRES solves A e = f − A u0 for the correction e from e = 0.
    P = SchwarzSweepPreconditioner(dd, alg.sweep, alg.coarse)
    e, stats = gmres(_RowOperator(dd.At), r0; N=P, memory=alg.memory, restart=true, atol=zero(T),
        rtol=T(alg.tol) * _rhs_scale(prob) / norm(r0), itmax=alg.maxiter, history=true)
    u = prob.u0 + e
    residuals = T.(stats.residuals) ./ _rhs_scale(prob)
    residuals[end] = _relative_residual(prob, u)
    return SchwarzResult(u, localize(dd, u), residuals, zeros(T, length(residuals)),
        stats.niter, last(residuals) <= alg.tol)
end

# ===== BiCGStab acceleration =====

"""
    SchwarzBiCGStab(; sweep=ParallelSweep(), coarse=nothing, tol=1e-8, maxiter=1000)

BiCGStab (van der Vorst 1992) right-preconditioned with one step of the
stationary Schwarz iteration ([`SchwarzSweepPreconditioner`](@ref)), for
nonsymmetric ``A``. Each iteration applies ``A`` and the preconditioner twice
and needs only four inner products, against GMRES's growing orthogonalization;
every vector operation is threaded, so with a threaded preconditioner no part
of an iteration runs on one core. Iteration stops when
``\\|f - A u\\| \\le \\mathrm{tol}\\,\\|f\\|``; the relative residual of the final
iterate is checked and recorded. A breakdown (vanishing ``\\hat r^\\mathsf{T} r``
or ``\\hat r^\\mathsf{T} v``) restarts from the current iterate.
"""
Base.@kwdef struct SchwarzBiCGStab{S<:SchwarzSweep,C<:Union{Nothing,AbstractCoarseSpace}}
    sweep::S = ParallelSweep()
    coarse::C = nothing
    tol::Float64 = 1e-8
    maxiter::Int = 1000
end

function solve(prob::SchwarzProblem{T}, alg::SchwarzBiCGStab) where {T}
    @argcheck alg.maxiter >= 1
    dd = prob.decomposition
    P = SchwarzSweepPreconditioner(dd, alg.sweep, alg.coarse)
    u = copy(prob.u0)
    iterations, _, residuals = _bicgstab!(u, (y, x) -> _threaded_mul!(y, dd.At, x), (y, x) -> mul!(y, P, x),
        prob.rhs; tol=alg.tol, maxiter=alg.maxiter)
    residuals ./= _rhs_scale(prob)
    residuals[end] = _relative_residual(prob, u)
    return SchwarzResult(u, localize(dd, u), residuals, zeros(T, length(residuals)),
        iterations, last(residuals) <= alg.tol)
end

# Contiguous index ranges for threaded vector operations.
function _chunks(n::Int)
    k = clamp(cld(n, 8192), 1, 4 * Threads.nthreads())
    return [((c - 1) * n ÷ k + 1):(c * n ÷ k) for c in 1:k]
end

# Threaded Σ f(i) over 1:n, summed in a fixed order (deterministic).
function _threaded_sum(f, n::Int, ::Type{T}) where {T}
    ranges = _chunks(n)
    partial = zeros(T, length(ranges))
    Threads.@threads for c in eachindex(ranges)
        acc = zero(T)
        @inbounds @simd for i in ranges[c]
            acc += f(i)
        end
        partial[c] = acc
    end
    return sum(partial)
end

function _threaded_foreach(f, n::Int)
    ranges = _chunks(n)
    Threads.@threads for c in eachindex(ranges)
        @inbounds for i in ranges[c]
            f(i)
        end
    end
end

"""
    _bicgstab!(x, A!, M!, b; tol, maxiter) -> (iterations, converged, residual norms)

Right-preconditioned BiCGStab for A x = b from the initial `x`, overwritten
with the solution. `A!(y, v)` sets y = A v and `M!(y, v)` sets y = M⁻¹ v.
Stops when ‖b − A x‖ ≤ tol ‖b‖. All vector operations are threaded.
"""
function _bicgstab!(x::Vector{T}, A!, M!, b::AbstractVector; tol::Real, maxiter::Integer) where {T}
    n = length(b)
    r, rhat, p, v, s, t, phat, shat = (zeros(T, n) for _ in 1:8)
    A!(r, x)
    _threaded_foreach(i -> (r[i] = b[i] - r[i]), n)
    copyto!(rhat, r)
    target = tol * max(sqrt(_threaded_sum(i -> abs2(b[i]), n, T)), eps(T))
    residuals = T[sqrt(_threaded_sum(i -> abs2(r[i]), n, T))]
    last(residuals) <= target && return 0, true, residuals
    ρ = α = ω = one(T)
    for it in 1:maxiter
        ρnew = _threaded_sum(i -> rhat[i] * r[i], n, T)
        if iszero(ρnew)                               # breakdown: restart the shadow residual
            copyto!(rhat, r); fill!(p, zero(T)); fill!(v, zero(T))
            ρ = α = ω = one(T)
            ρnew = _threaded_sum(i -> abs2(r[i]), n, T)
        end
        β = (ρnew / ρ) * (α / ω)
        let β = β, ω = ω                           # (let: closures must not capture reassigned locals)
            _threaded_foreach(i -> (p[i] = r[i] + β * (p[i] - ω * v[i])), n)
        end
        M!(phat, p)
        A!(v, phat)
        σ = _threaded_sum(i -> rhat[i] * v[i], n, T)
        if iszero(σ)
            copyto!(rhat, r); fill!(p, zero(T)); fill!(v, zero(T))
            ρ = α = ω = one(T)
            push!(residuals, last(residuals))
            continue
        end
        α = ρnew / σ
        let α = α
            _threaded_foreach(i -> (s[i] = r[i] - α * v[i]), n)
        end
        snorm = sqrt(_threaded_sum(i -> abs2(s[i]), n, T))
        if snorm <= target
            let α = α
                _threaded_foreach(i -> (x[i] += α * phat[i]), n)
            end
            push!(residuals, snorm)
            return it, true, residuals
        end
        M!(shat, s)
        A!(t, shat)
        ω = _threaded_sum(i -> t[i] * s[i], n, T) / _threaded_sum(i -> abs2(t[i]), n, T)
        let α = α, ω = ω
            _threaded_foreach(i -> (x[i] += α * phat[i] + ω * shat[i]; r[i] = s[i] - ω * t[i]), n)
        end
        ρ = ρnew
        push!(residuals, sqrt(_threaded_sum(i -> abs2(r[i]), n, T)))
        last(residuals) <= target && return it, true, residuals
    end
    return maxiter, false, residuals
end
