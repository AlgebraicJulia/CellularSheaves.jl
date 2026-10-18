# Sheaf ADMM for SchwarzMethods: the Schwarz cover as a consensus problem over
# its overlap sheaf. Included in SchwarzMethods.

export SheafADMM, LocalObjective, local_objectives, ADMMPenalty, StalkPenalty, SharedPenalty

# ===== Sheaf ADMM =====

"""
    LocalObjective

The quadratic local objective ``f_i(x_i) = \\tfrac12 x_i^\\mathsf{T} K_i x_i - b_i^\\mathsf{T} x_i``
of one subdomain, on its vertex stalk ``\\mathbb R^{\\overline\\Omega_i}``, with
`matrix` ``K_i`` and `rhs` ``b_i``. [`local_objectives`](@ref) builds them so
that on every global section ``x = \\mathrm{localize}(u)``

```math
\\sum_i f_i(x_i) = \\tfrac12\\, u^\\mathsf{T} A u - f^\\mathsf{T} u .
```
"""
struct LocalObjective{T}
    matrix::SparseMatrixCSC{T,Int}
    rhs::Vector{T}
end

"""
    local_objectives(prob::SchwarzProblem) -> Vector{LocalObjective}

Split the energy ``\\tfrac12 u^\\mathsf{T} A u - f^\\mathsf{T} u`` of `prob` into one
[`LocalObjective`](@ref) per vertex stalk of the ghost-layer cover. The energy is
a sum of edge and vertex terms,

```math
u^\\mathsf{T} A u = \\sum_{k<l} |a_{kl}|\\,(u_k + \\operatorname{sign}(a_{kl})\\, u_l)^2
    + \\sum_k s_k\\, u_k^2, \\qquad s_k = a_{kk} - \\sum_{l \\neq k} |a_{kl}| .
```

Each edge term is shared equally among the stalks containing both ``k`` and
``l``, and each vertex term and ``f_k`` among the stalks containing ``k``. Every
edge of `A` lies in some closed subdomain, so nothing is lost. When `A` is
diagonally dominant (``s_k \\ge 0``, e.g. an M-matrix) every term is a square,
so every ``K_i`` is positive semidefinite.
"""
function local_objectives(prob::SchwarzProblem{T}) where {T}
    dd = prob.decomposition
    S, cover = dd.A, dd.cover
    entries = [(Int[], Int[], T[]) for _ in cover.stalks]
    rhs = [zeros(T, length(s)) for s in cover.stalks]
    function add!(i, a, b, v)
        push!(entries[i][1], a); push!(entries[i][2], b); push!(entries[i][3], v)
    end
    for k in axes(S, 2)
        holders = cover.members[k]
        excess = S[k, k] - sum(abs(nonzeros(S)[p]) for p in nzrange(S, k) if rowvals(S)[p] != k; init=zero(T))
        for (i, ℓ) in holders
            add!(i, ℓ, ℓ, excess / length(holders))
            rhs[i][ℓ] += prob.rhs[k] / length(holders)
        end
        for p in nzrange(S, k)
            l, a = rowvals(S)[p], nonzeros(S)[p]
            (l > k && !iszero(a)) || continue
            common = [(i, ℓk, ℓl) for (i, ℓk) in holders for (j, ℓl) in cover.members[l] if i == j]
            w, σ = abs(a) / length(common), sign(a)
            for (i, ℓk, ℓl) in common
                add!(i, ℓk, ℓk, w); add!(i, ℓl, ℓl, w)
                add!(i, ℓk, ℓl, σ * w); add!(i, ℓl, ℓk, σ * w)
            end
        end
    end
    return [LocalObjective(sparse(I, J, V, length(s), length(s)), b)
            for ((I, J, V), b, s) in zip(entries, rhs, cover.stalks)]
end

"""
    ADMMPenalty

Which dofs of each stalk the penalty ``\\tfrac{\\rho}{2}\\lVert x_i - z_i + y_i\\rVert_{P_i}^2``
of [`SheafADMM`](@ref) acts on: [`StalkPenalty`](@ref) or [`SharedPenalty`](@ref).
"""
abstract type ADMMPenalty end

"""
    StalkPenalty()

``P_i = I``: the penalty acts on the whole stalk (Hanks et al. 2025, Algorithm 1).
"""
struct StalkPenalty <: ADMMPenalty end

"""
    SharedPenalty()

``P_i`` selects the dofs held by more than one stalk (the overlaps and ghost
layers), where the consensus constraint binds.
"""
struct SharedPenalty <: ADMMPenalty end

_penalty_mask(::StalkPenalty, cover, stalk) = trues(length(stalk))
_penalty_mask(::SharedPenalty, cover, stalk) = BitVector(length(cover.members[k]) > 1 for k in stalk)

"""
    SheafADMM(; rho, penalty=StalkPenalty(), projection_steps=nothing, tol=1e-8, maxiter=1000)

Sheaf ADMM (Hanks, Riess, Cohen, Gross, Hale, Fairbanks 2025, Algorithm 1) on
the ghost-layer overlap sheaf. It solves the homological program

```math
\\min_{x \\in C^0} \\sum_i f_i(x_i) \\quad \\text{subject to} \\quad x \\in H^0 ,
```

with the [`local_objectives`](@ref) of the PDE, so its solution glues to the
solution of ``A u = f``. With copies ``z`` and scaled multipliers ``y`` on the
vertex stalks, each iteration is

```math
\\begin{aligned}
x_i &\\leftarrow \\operatorname{argmin}_{x_i} f_i(x_i) + \\tfrac{\\rho}{2}\\lVert x_i - z_i + y_i\\rVert_{P_i}^2
     = (K_i + \\rho P_i)^{-1}\\big(b_i + \\rho P_i (z_i - y_i)\\big),\\\\
z &\\leftarrow \\Pi_{H^0}(x + y),\\\\
y_i &\\leftarrow y_i + x_i - z_i .
\\end{aligned}
```

- `penalty = StalkPenalty()` is the method of the paper: ``P_i = I``, the penalty acts on
  the whole stalk. [`SharedPenalty`](@ref)`()` penalizes only the dofs held by more than
  one stalk (the overlaps and ghost layers), which is where the consensus
  constraint actually binds. It is closer to a Robin transmission condition,
  whose penalty acts only on the interface (see the comparison in the feature
  guide).
- ``\\Pi_{H^0}`` is the orthogonal projection onto global sections. For this
  sheaf it averages the copies of every dof, one exchange with the neighbours
  (`projection_steps = nothing`). With `projection_steps = k` it is replaced by
  ``k`` explicit sheaf-diffusion steps ``z \\leftarrow z - \\alpha L_F z`` with
  ``\\alpha = 1/c_{\\max}`` (``c_{\\max}`` the largest number of copies of a dof),
  as in Algorithm 1, at one neighbour exchange per step.

The local matrices ``K_i + \\rho P_i`` are factored once with `ChordalLDLt`.
Convergence (Boyd et al. 2011, §3.2) needs convex ``f_i`` and an exact
projection, so `A` should be diagonally dominant; the local matrices are
checked for positive definiteness.
"""
Base.@kwdef struct SheafADMM{P<:ADMMPenalty}
    rho::Float64
    penalty::P = StalkPenalty()
    projection_steps::Union{Nothing,Int} = nothing
    tol::Float64 = 1e-8
    maxiter::Int = 1000
end

# Π_{H^0}: average the copies of every dof.
function _project_sections!(zs, dd::SchwarzDecomposition{T}, vs) where {T}
    total = zeros(T, size(dd.A, 1))
    copies = zeros(Int, size(dd.A, 1))
    for (v, s) in zip(vs, dd.cover.stalks)
        total[s] .+= v
        copies[s] .+= 1
    end
    total ./= copies
    for (z, s) in zip(zs, dd.cover.stalks)
        z .= view(total, s)
    end
    return zs
end

# k explicit sheaf-diffusion steps z ← z − α L z on the overlap sheaf.
function _diffuse_sections!(zs, dd::SchwarzDecomposition, vs, steps::Int)
    cover = dd.cover
    α = 1 / maximum(length, cover.members)
    foreach(copyto!, zs, vs)
    for _ in 1:steps
        Lz = [zero(z) for z in zs]
        for ((i, j), mine) in cover.overlaps
            theirs = cover.overlaps[j => i]
            Lz[i][mine] .+= zs[i][mine] .- zs[j][theirs]
        end
        foreach((z, l) -> (z .-= α .* l), zs, Lz)
    end
    return zs
end

function solve(prob::SchwarzProblem{T}, alg::SheafADMM) where {T}
    @argcheck issymmetric(prob.decomposition.A) "sheaf ADMM needs a symmetric A"
    @argcheck alg.rho > 0 "the ADMM penalty must be positive"
    @argcheck alg.projection_steps === nothing || alg.projection_steps >= 1
    dd = prob.decomposition
    cover = dd.cover
    objectives = local_objectives(prob)
    masks = [_penalty_mask(alg.penalty, cover, s) for s in cover.stalks]
    factors = map(enumerate(objectives)) do (i, obj)
        Ki = obj.matrix + alg.rho * spdiagm(T.(masks[i]))
        F = ldlt!(ChordalLDLt(Ki), RowMaximum(); check=false)
        @argcheck all(>(0), F.D.diag) "the ADMM local matrix of subdomain $i is not positive definite; increase rho or use penalty=StalkPenalty()"
        F
    end

    z = localize(dd, prob.u0)
    zs = blocks(z)
    xs = [copy(zi) for zi in zs]
    ys = [zero(zi) for zi in zs]
    u = glue(dd, z)
    residuals = T[_relative_residual(prob, u)]
    disagreements = T[0]
    iterations = 0
    while last(residuals) > alg.tol && iterations < alg.maxiter
        Threads.@threads for i in eachindex(xs)
            b = objectives[i].rhs .+ alg.rho .* masks[i] .* (zs[i] .- ys[i])
            xs[i] .= _ldlt_solve(factors[i], b)
        end
        vs = [x .+ y for (x, y) in zip(xs, ys)]
        alg.projection_steps === nothing ? _project_sections!(zs, dd, vs) :
            _diffuse_sections!(zs, dd, vs, alg.projection_steps)
        for (y, x, zi) in zip(ys, xs, zs)
            y .+= x .- zi
        end
        iterations += 1
        u = glue(dd, z)
        push!(residuals, _relative_residual(prob, u))
        push!(disagreements, overlap_disagreement(dd, mortar(xs)))
    end
    return SchwarzResult(u, mortar(xs), residuals, disagreements, iterations, last(residuals) <= alg.tol)
end
