@kwdef struct HSDSettings{T} <: AbstractSettings{T}
    verbose::Int = 0
    step_frac::T = 0.99
    feas_tol::T = 1e-8
    gap_tol::T = 1e-8
    itmax::Int = 100
    near_factor::T = 1000.0
    stall_tol::T = 1e-6
    forcing_frac::T = 0.1
    forcing_ceil::T = 0.3
    refine_itmax::Int = 3
    refine_stall::T = 0.5
    floor_patience::Int = 3
    scale_itmax::Int = 10
    rgmin::T = 1e-9
    rgmax::T = 1e-6
    kkt::KKTSettings{T} = UzawaSettings{T}()
    # HSD-specific
    illposed_tol::T = 1e-10
    infeas_abs::T = 1e-8
    infeas_rel::T = 1e-8
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
    @printf(io, "%srgmin:        %8.2e  rgmax:        %8.2e\n", pad, set.rgmin, set.rgmax)
    @printf(io, "%sillposed_tol: %8.2e  infeas_abs:   %8.2e\n", pad, set.illposed_tol, set.infeas_abs)
    @printf(io, "%sinfeas_rel:   %8.2e\n", pad, set.infeas_rel)
    println(io)
    println(io, pad, "kkt: ", typeof(set.kkt))
    showsettings(io, set.kkt; indent=indent+4)
    return
end

