@kwdef struct HSDSettings{T, E <: EliminationAlgorithm} <: AbstractSettings{T}
    verbose::Int = 0
    step_frac::T = 0.99
    feas_tol::T = 1e-8
    gap_tol::T = 1e-8
    max_iter::Int = 100
    near_factor::T = 1000.0
    stall_tol::T = 1e-3         # min per-iteration relative progress; stall = below this for a window
    refine_stall_tol::T = 0.5   # min per-pass residual reduction in the KKT refinement loop
    scale_max_iter::Int = 10    # equilibration sweeps
    refine_max_iter::Int = 10   # max KKT iterative-refinement passes
    newton_max_iter::Int = 100  # max CG iterations per KKT (Newton-system) solve
    aug_tol::T = 1e-7           # relative augmentation, in initial-problem units (the anchored knob)
    pivot::Bool = false         # rank-revealing pivoted base factorization
    elim_alg::E = DEFAULT_ELIMINATION_ALGORITHM   # KKT elimination ordering (construction-time)
    # HSD-specific
    illposed_tol::T = 1e-10
    infeas_abs::T = 1e-8
    infeas_rel::T = 1e-8
end

# infer E from the elim_alg argument so HSDSettings{T}(...) stays the call form
function HSDSettings{T}(; elim_alg::E = DEFAULT_ELIMINATION_ALGORITHM, kwargs...) where {T, E <: EliminationAlgorithm}
    return HSDSettings{T, E}(; elim_alg, kwargs...)
end

function HSDSettings{T}(set::HSDSettings; kw...) where {T}
    nt = (fn => getfield(set, fn) for fn in fieldnames(HSDSettings))
    return HSDSettings{T}(; nt..., kw...)
end

function HSDSettings(set::HSDSettings{T}; kw...) where {T}
    return HSDSettings{T}(set; kw...)
end

function showsettings(io::IO, set::HSDSettings; indent::Integer=0)
    pad = " "^indent
    @printf(io, "%sverbose:      %8s  feas_tol:      %8.2e\n", pad, set.verbose, set.feas_tol)
    @printf(io, "%sstep_frac:    %8.2e  gap_tol:       %8.2e\n", pad, set.step_frac, set.gap_tol)
    @printf(io, "%smax_iter:     %8d  stall_tol:     %8.2e\n", pad, set.max_iter, set.stall_tol)
    @printf(io, "%snear_factor:  %8.2e  refine_stall_tol: %5.2e\n", pad, set.near_factor, set.refine_stall_tol)
    @printf(io, "%sscale_max_iter:%7d  refine_max_iter: %6d\n", pad, set.scale_max_iter, set.refine_max_iter)
    @printf(io, "%snewton_max_iter:%6d  aug_tol:       %8.2e\n", pad, set.newton_max_iter, set.aug_tol)
    @printf(io, "%sillposed_tol: %8.2e  infeas_abs:    %8.2e\n", pad, set.illposed_tol, set.infeas_abs)
    @printf(io, "%sinfeas_rel:   %8.2e\n", pad, set.infeas_rel)
    return
end
