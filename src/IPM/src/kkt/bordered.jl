# Bordered KKT solver — wraps a 2-row KKTSolver and adds the HSD τ-border via a Woodbury/Schur update.
# The Woodbury column w₂ = KKT⁻¹[c; g] and the capacitance scalar S are cached: they are fixed for a
# given factorization and reused by the predictor and corrector solves at that factorization. The border
# data (c, g, Q, Hc) is passed per call, not stored. initkkt! factors + solves/caches w₂ + capacitance;
# solvekkt! does the bordered base solve (newton!) + the 3-row refinement (refinehsd!).
struct BorderedSolver{T, W <: KKTSolver{T}} <: KKTSolver{T}
    inner::W          # the 2-row solver (factor + base solve primitive)
    Δp2::Vector{T}    # cached Woodbury column, primal part (n)
    Δy2::Vector{T}    # cached Woodbury column, dual part (m)
    aτ::Vector{T}     # cached border row c - 2Qp/τ (n)
    QΔp2::Vector{T}   # capacitance scratch (n)
    S::Scalar{T}      # cached capacitance scalar
end

function BorderedSolver(inner::W, B::BlockSparseMatrix{T, I}) where {T, I, W <: KKTSolver{T}}
    m, n = size(B)
    return BorderedSolver(inner, zeros(T, n), zeros(T, m), zeros(T, n), zeros(T, n), ones(T))
end
