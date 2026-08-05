@kwdef struct HSDSettings{T, E <: EliminationAlgorithm} <: AbstractSettings{T}
    verbose::Int = 0
    step_frac::T = 0.99
    feas_tol::T = 1e-8
    gap_tol::T = 1e-8
    itmax::Int = 100
    near_factor::T = 1000.0
    stall_tol::T = 1e-6
    forcing_frac::T = 0.1
    forcing_ceil::T = 0.3
    refine_itmax::Int = 10   # final_spec.md C1: raised from 3 — un-gated refinement needs the budget
    refine_stall::T = 0.5
    floor_patience::Int = 3
    scale_itmax::Int = 10
    aaug::T = zero(T)   # absolute augmentation (raw-α override: set raug=0, aaug=α_raw)
    raug::T = 1e7       # relative augmentation, in initial-problem units (the anchored knob)
    fix_alpha::Bool = false   # ORACLE ONLY: setaug! becomes a no-op so a caller-injected α survives step!
    policy::Int = 0     # α policy: 0 = current (fudged window + cap), 1 = policies_1_and_2.md Policy 1
    vartol::Bool = false   # DEPRECATED (maps to forcing=1 when forcing=0); use `forcing`
    forcing::Int = 0    # 0=absolute μ-schedule, 1=vartol (η·R0), 2=gnorm (δ·‖g‖), 4=fixtol (floor)
    tol0::T = 1e-3      # vartol: relative reduction target (Z&G tune 1e-3..1e-2 per problem class)
    delta::T = 0.05     # gnorm: relative primal-residual reduction ‖g−BΔp‖ ≤ δ·‖g‖ (fixed constant)
    abstol::T = 1e-8    # fixtol: the fixed (scaled) residual target solved to every iteration
    elim::E = DEFAULT_ELIMINATION_ALGORITHM   # KKT elimination ordering (construction-time)
    # HSD-specific
    illposed_tol::T = 1e-10
    infeas_abs::T = 1e-8
    infeas_rel::T = 1e-8
end

# infer E from the elim argument so HSDSettings{T}(...) stays the call form
function HSDSettings{T}(; elim::E = DEFAULT_ELIMINATION_ALGORITHM, kwargs...) where {T, E <: EliminationAlgorithm}
    return HSDSettings{T, E}(; elim, kwargs...)
end

function showsettings(io::IO, set::HSDSettings; indent::Integer=0)
    pad = " "^indent
    @printf(io, "%sverbose:      %8s  feas_tol:     %8.2e\n", pad, set.verbose, set.feas_tol)
    @printf(io, "%sstep_frac:    %8.2e  gap_tol:      %8.2e\n", pad, set.step_frac, set.gap_tol)
    @printf(io, "%sitmax:        %8d  stall_tol:    %8.2e\n", pad, set.itmax, set.stall_tol)
    @printf(io, "%snear_factor:  %8.2e  forcing_frac: %8.2e\n", pad, set.near_factor, set.forcing_frac)
    @printf(io, "%sforcing_ceil: %8.2e\n", pad, set.forcing_ceil)
    @printf(io, "%srefine_itmax: %8d  refine_stall: %8.2e\n", pad, set.refine_itmax, set.refine_stall)
    @printf(io, "%sscale_itmax:  %8d\n", pad, set.scale_itmax)
    @printf(io, "%sillposed_tol: %8.2e  infeas_abs:   %8.2e\n", pad, set.illposed_tol, set.infeas_abs)
    @printf(io, "%sinfeas_rel:   %8.2e\n", pad, set.infeas_rel)
    @printf(io, "%saaug:         %8.2e  raug:         %8.2e\n", pad, set.aaug, set.raug)
    return
end

