# =============================================================================
# oracle.jl — parsimonious per-iteration α oracle (measurement infra; deleted at end of campaign)
# =============================================================================
#
# At each IPM/HSD iteration, from the current solver state, take one REAL production step at every
# candidate α (deepcopy the solver, set α, step!), score it, advance along the best α. Each candidate
# is a genuine warm-started step on a full deepcopy — production truth, same warm start / factorization
# / two-or-three-solve structure the live solver pays.
#
# Score = (state, ncraig) compared lexicographically:
#   • state  = 0 iff EVERY role's refinement reached force/floor this step, else 1 (any itmax/stall)
#   • ncraig = summed base+refn CRAIG over the step's solves (woodbury role included on HSD)
# Ties resolve to the LOWEST α (grid ascending, strict <).
#
# Records one lean NamedTuple per (iter, α): the cost object (ncraig = Σ_roles pbase+prefn, which is
# exactly the theory's (1+j)+n0+Σni per solve, summed over the step's solves), its decomposition
# (per-role base / refn CRAIG + pass counts + statuses), the two FREE scales' raw inputs
# (pfres=‖r0‖ primal base residual → ℓ_p=α‖r0‖; pdres≈s0 dual → ℓ_d=s0/α), ρ, μ, tolerances, ‖B‖.
# No spectral harvest, no per-pass trace — those are added only if/when the policy needs them (W3).

using LinearAlgebra: norm, Symmetric
using Printf: @sprintf
using CellularSheaves.IPM: step!, mu, CONTINUE, REACHED_FORCE, REACHED_FLOOR,
    getaug, CTRL_CAP_IPM, CTRL_CAP_HSD, IPMSolver, HSDSolver

# The controller cap the live setaug! would use, per solver kind.
_ctrl_cap(::IPMSolver) = CTRL_CAP_IPM
_ctrl_cap(::HSDSolver) = CTRL_CAP_HSD

const DEFAULT_ALPHA_GRID = [round(10.0^e, sigdigits=4) for e in 0.0:0.5:18.0]   # half-decades 1e0..1e18 (37 pts)

_reached(st) = st === REACHED_FORCE || st === REACHED_FLOOR

function _oracle_state(row)
    r = _reached(row.pstat) && _reached(row.cstat)
    hasproperty(row, :wstat) && (r &= _reached(row.wstat))
    return r ? 0 : 1
end

function _oracle_craig(row)
    c = row.pbase + row.prefn + row.cbase + row.crefn
    hasproperty(row, :wbase) && (c += row.wbase + row.wrefn)
    return c
end

function _oracle_record(i, α, sc, st, αctrl)
    row = sc.hist[end]
    μ   = row.μ
    μ1  = first(sc.hist.μ)
    T   = eltype(sc.wrk.rp)
    # exact recomputation of the tolerances step! used this iteration
    ftol  = min(sc.settings.forcing_frac * μ / μ1, sc.settings.forcing_ceil)
    fltol = 100 * eps(T) * (one(T) + max(norm(sc.wrk.rp, Inf), norm(sc.wrk.rd, Inf)))
    hasw  = hasproperty(row, :wstat)
    return (iter = i, alpha = α, alpha_ctrl = αctrl, state = _oracle_state(row), ncraig = _oracle_craig(row),
            ipm_status = st, mu = μ, mu_next = mu(sc), rho = row.ρ,
            force_tol = ftol, floor_tol = fltol, normB = norm(sc.B),
            step = row.step, ipm_pres = row.pres, ipm_dres = row.dres,
            # cost decomposition, per role (base includes the "1" backsolve; refn = j+Σni)
            pbase = row.pbase, prefn = row.prefn, ppass = row.ppass, pstat = row.pstat,
            cbase = row.cbase, crefn = row.crefn, cpass = row.cpass, cstat = row.cstat,
            wbase = hasw ? row.wbase : nothing, wrefn = hasw ? row.wrefn : nothing,
            wpass = hasw ? row.wpass : nothing, wstat = hasw ? row.wstat : nothing,
            # the two FREE scales' raw inputs: base residual (‖r0‖) + pass-1 entry residuals (≈ s0)
            pfres = row.pfres, ppres = row.ppres, pdres = row.pdres,
            cfres = row.cfres, cpres = row.cpres, cdres = row.cdres,
            wfres = hasw ? row.wfres : nothing, wpres = hasw ? row.wpres : nothing,
            wdres = hasw ? row.wdres : nothing,
            chosen = false)
end

# Serialize records to CSV. Enums stringified, `nothing`→empty, floats scientific (NaN preserved).
_csvcell(::Nothing) = ""
_csvcell(x::AbstractFloat) = isnan(x) ? "NaN" : @sprintf("%.15e", x)
_csvcell(x) = string(x)
function write_oracle_csv(path::AbstractString, records::AbstractVector{<:NamedTuple})
    open(path, "w") do io
        isempty(records) && return
        cols = keys(records[1])
        println(io, join(cols, ","))
        for r in records
            println(io, join((_csvcell(getproperty(r, c)) for c in cols), ","))
        end
    end
    return path
end

function solve_logged(s0, grid = DEFAULT_ALPHA_GRID; itmax::Integer = s0.settings.itmax)
    s = s0
    records = NamedTuple[]
    gs = sort(collect(grid))                       # ascending ⇒ ties on (state,ncraig) → lowest α
    i = 0
    while i < itmax
        i += 1
        iter_recs = NamedTuple[]
        best = nothing
        bestscore = (typemax(Int), typemax(Int))
        beststatus = CONTINUE
        bestidx = 0
        # α the live controller (setaug!) would have picked at THIS iterate — the oracle disables the
        # controller (fix_alpha) so the swept α survives, but we still record its choice to judge it.
        # iter≥2: exact — getaug reads only the previous step's recorded α/αmin/αmax. iter 1: raug
        # formula, filled below from a stepped candidate's scaled H (scaling is α-independent).
        T = eltype(s.wrk.rp)
        αctrl = isempty(s.hist) ? T(NaN) : getaug(s.hist, _ctrl_cap(s))
        for α in gs
            sc = deepcopy(s)
            sc.α[] = α
            st = step!(sc)
            isnan(αctrl) && (αctrl = s.settings.aaug + s.settings.raug * norm(Symmetric(sc.H, :L)) / s.nB[]^2)
            push!(iter_recs, _oracle_record(i, α, sc, st, αctrl))
            score = (iter_recs[end].state, iter_recs[end].ncraig)
            if score < bestscore                   # strict: first (lowest) α at the best score wins
                bestscore = score
                best = sc
                beststatus = st
                bestidx = lastindex(iter_recs)
            end
        end
        # entry already converged ⇒ every candidate trivially returns ncraig=0 at a terminal status;
        # record nothing and stop (removes the artifact).
        if beststatus !== CONTINUE && iter_recs[bestidx].ncraig == 0
            break
        end
        iter_recs[bestidx] = merge(iter_recs[bestidx], (chosen = true,))
        append!(records, iter_recs)
        s = best
        beststatus === CONTINUE || break           # terminal status at the chosen α ⇒ done
    end
    return s, records
end
