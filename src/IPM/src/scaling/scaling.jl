#
# Static Ruiz equilibration (preprocessing)
# ==========================================
#
# Equilibrates the conic QP
#
#   primal: min f'p + ½ p'Q p   s.t.  B p = g,  p ∈ K
#   dual:   max …                s.t.  B'y + d − Q p = f,  d ∈ K*
#
# by a block-constant column scaling D and a per-row scaling E so that, after
# scaling, every constraint row of B and every cone-block column of B has
# ∞-norm ≈ 1. This is the classic Ruiz (2001) row/column equilibration,
# specialized so that it is compatible with conic geometry.
#
# The scaled data is
#
#   B̂ = E B D,   ĝ = E g,   ĉ = D f,   Q̂ = D Q D
#
# The change of variables is
#
#   p̂ = D⁻¹ p,   d̂ = D d,   ŷ = E⁻¹ y
#
# so that a scaled solution (p̂, d̂, ŷ) maps back to the original problem by
#
#   p = D p̂,   d = D⁻¹ d̂,   y = E ŷ.
#
# Why one scalar per cone block.
# ------------------------------
# D must be CONSTANT on each cone block: D = blkdiag(t_v · I) with t_v > 0. The
# columns of a block are the svec/coordinates of a single conic variable, and
# every cone in this package (POS, SOC, SDP, EXP, POW, and the free cone) is
# closed under positive *uniform* scaling — t·x ∈ K ⇔ x ∈ K — but NOT under
# per-coordinate scaling (e.g. shrinking only x₀ of a second-order cone, or only
# the off-diagonal svec entries of a PSD block, leaves the cone). A single
# positive scalar per block is therefore the most aggressive column scaling that
# is guaranteed to satisfy  p ∈ K ⇔ D⁻¹ p ∈ K  for *all* cone types. Rows carry
# no cone (the dual y is free), so E is unconstrained per row.
#
# The column scaling acts symmetrically on Q as Q̂ = D Q D with D constant on each
# block, so an off-diagonal block Q[u,v] is scaled by vscl[u]·vscl[v].
#

abstract type AbstractScaling{T} end

include("ipm.jl")
include("hsd.jl")

#
# Compute row and vertex-block ∞-norms in a single pass.
# rnrm[i] ← max over row i of B
# vnrm[v] ← max over vertex block v of B and Q
#
function infnorms!(rnrm::AbstractVector, vnrm::AbstractVector, B::BlockSparseMatrix, Q::BlockSparseMatrix)
    fill!(rnrm, false)

    for v in vtxs(B)
        nv = zero(eltype(vnrm))

        for e in srcrange(Q, v)
            nv = max(nv, norm(block(Q, Q.tgt[e], v, e), Inf))
        end

        for e in srcrange(B, v)
            u = B.tgt[e]
            ru = rowrange(B, u)
            Be = block(B, u, v, e)

            for jloc in axes(Be, 2)
               for iloc in axes(Be, 1)
                    i   = ru[iloc]
                    Bij = Be[iloc, jloc]

                    nv      = max(nv,      abs(Bij))
                    rnrm[i] = max(rnrm[i], abs(Bij))
                end
            end
        end

        vnrm[v] = nv
    end

    return
end

#
# B ← E B D  and  Q ← D Q D, where E = diag(rscl) (per row) and D = diag with the
# block-constant value vscl[v] on block v. Applied in place to the working copies.
#
function applyscaling!(B::BlockSparseMatrix, Q::BlockSparseMatrix, rscl::AbstractVector, vscl::AbstractVector)
    for v in vtxs(B)
        sv = vscl[v]

        for e in srcrange(B, v)
            u  = B.tgt[e]
            ru = rowrange(B, u)
            Be = block(B, u, v, e)

            for jloc in axes(Be, 2)
                for iloc in axes(Be, 1)
                    Be[iloc, jloc] *= rscl[ru[iloc]] * sv
                end
            end
        end

        for e in srcrange(Q, v)
            lmul!(vscl[Q.tgt[e]] * sv, block(Q, Q.tgt[e], v, e))
        end
    end

    return
end

"""
    equilibrate!(scaling, B, Q, f, g; itmax=10, tol=1e-3)

Run block-aware Ruiz equilibration on the conic-QP data in place.
Mutates `scaling`, `B`, `Q`, `f`, and `g`. Returns `scaling`.

Keyword arguments:
- `itmax` : maximum Ruiz sweeps.
- `tol`   : stop when every row/block ∞-norm is within `tol` of 1.
"""
function equilibrate!(
        scaling::AbstractScaling{T},
        B::BlockSparseMatrix{T},
        Q::BlockSparseMatrix{T},
        f::AbstractVector{T},
        g::AbstractVector{T};
        itmax::Int=10,
        tol::Real=1e-3
    ) where {T}
    nrow = nrows(B)
    nvtx = nvtxs(B)

    cscl = scaling.cscl
    rscl = scaling.rscl

    vscl = ones(T, nvtx)
    vnrm = zeros(T, nvtx)
    rnrm = zeros(T, nrow)
    vswp = ones(T, nvtx)
    rswp = ones(T, nrow)

    for _ in 1:itmax
        infnorms!(rnrm, vnrm, B, Q)

        d = zero(T)

        for v in vtxs(B)
            n = sqrt(vnrm[v])
            d = max(d, abs(one(T) - n))

            if !iszero(n)
                vscl[v] *= vswp[v] = inv(n)
            end
        end

        for i in rows(B)
            n = sqrt(rnrm[i])
            d = max(d, abs(one(T) - n))

            if !iszero(n)
                rscl[i] *= rswp[i] = inv(n)
            end
        end

        applyscaling!(B, Q, rswp, vswp)

        d < tol && break
    end

    for v in vtxs(B)
        s = vscl[v]

        for j in colrange(B, v)
            f[j] *= cscl[j] = s
        end
    end

    for i in rows(B)
        g[i] *= rscl[i]
    end

    return scaling
end

"""
    scale!(p, d, y, scaling)

Apply forward scaling to vectors in place:

    p ← D⁻¹ p,   d ← D d,   y ← E⁻¹ y
"""
function scale!(p::AbstractVector, d::AbstractVector, y::AbstractVector, scaling::AbstractScaling)
    p ./= scaling.cscl
    d .*= scaling.cscl
    y ./= scaling.rscl
    return p, d, y
end

"""
    unscale!(p, d, y, scaling)

Apply inverse scaling to vectors in place:

    p ← D p,   d ← D⁻¹ d,   y ← E y
"""
function unscale!(p::AbstractVector, d::AbstractVector, y::AbstractVector, scaling::AbstractScaling)
    p .*= scaling.cscl
    d ./= scaling.cscl
    y .*= scaling.rscl
    return p, d, y
end
