############################################################################################
# solve_logged — greedy per-iteration α oracle (production units, no replay)
############################################################################################
#
# At each IPM iteration, from the current solver state, take one REAL production step at every
# candidate α (deepcopy the solver, set α, call step!), score it, and advance along the best α.
# Because each candidate is a genuine warm-started step on a full deepcopy, the cost is production
# truth — the same warm start / factorization / three-solve structure the live solver pays.
#
# Score is the tuple (state, ncraig), compared lexicographically:
#   • state  = 0 iff EVERY refinement status this step reached force/floor, else 1 (any itmax/stall)
#   • ncraig = summed base+refn CRAIG over the step's solves (woodbury role included on HSD)
# Ties resolve to the LOWEST α (grid is sorted ascending and the update uses a strict <).
#
# Returns (final_solver, records), where records is one NamedTuple per (iteration, α) carrying the
# FULL score — every per-role refinement status and CRAIG count — plus `chosen` marking the winner.

const DEFAULT_ALPHA_GRID = [round(10.0^e, sigdigits=4) for e in 0.0:0.5:18.0]   # half-decades 1e0..1e18 (37 pts)

# NOTE: σ̂²min/σ̂²max are no longer computed here. They are harvested INSIDE the predictor's base solve
# (`solve_uzw!` in kkt/uzawa.jl, via `gk_spectral`), on the exact rhs `r` that `craig!` consumes, and
# ride the history row as `s2min_p`/`s2max_p` — a stand-in for a future inline-CRAIG spectral harvest.

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

function _oracle_record(i, α, sc, st)
    row = sc.hist[end]
    μ   = row.μ
    μ1  = first(sc.hist.μ)
    T   = eltype(sc.wrk.rp)
    # exact recomputation of the tolerances step! used this iteration (μ + settings; persisted rp/rd)
    ftol  = min(sc.settings.forcing_frac * μ / μ1, sc.settings.forcing_ceil)
    fltol = 100 * eps(T) * (one(T) + max(norm(sc.wrk.rp, Inf), norm(sc.wrk.rd, Inf)))
    hasw  = hasproperty(row, :wstat)
    dg    = diag(sparse(sc.H))                     # scaled-Hessian diagonal (α-independent within an iter)
    Ld    = diag(sc.kkt.F)                          # Cholesky factor diagonal — κ(F) proxy / ρ-ladder early warning
    base = (iter = i, alpha = α, state = _oracle_state(row), ncraig = _oracle_craig(row),
            ipm_status = st, mu = μ, mu_next = mu(sc), rho = row.ρ,
            force_tol = ftol, floor_tol = fltol,
            sigma2min = row.s2min_p, sigma2max = row.s2max_p,   # harvested in the predictor base solve (rhs r)
            step = row.step,
            tau = hasproperty(row, :τ) ? row.τ : nothing,
            kappa = hasproperty(row, :κ) ? row.κ : nothing,
            norm_dp = norm(sc.wrk.Δp), norm_dy = norm(sc.wrk.Δy), norm_y = norm(sc.y),
            hdiag_min = minimum(dg), hdiag_max = maximum(dg),
            pstat = row.pstat, cstat = row.cstat, wstat = hasw ? row.wstat : nothing,
            pbase = row.pbase, prefn = row.prefn, ppass = row.ppass,
            cbase = row.cbase, crefn = row.crefn, cpass = row.cpass,
            wbase = hasw ? row.wbase : nothing, wrefn = hasw ? row.wrefn : nothing,
            wpass = hasw ? row.wpass : nothing,
            pres0 = row.pres0, pres1 = row.pres1, cres0 = row.cres0, cres1 = row.cres1,
            wres0 = hasw ? row.wres0 : nothing, wres1 = hasw ? row.wres1 : nothing,
            r0_p = row.r0_p, r0_c = row.r0_c, r0_w = hasw ? row.r0_w : nothing,
            craig_p = row.craig_p, craig_c = row.craig_c, craig_w = hasw ? row.craig_w : nothing,
            r1_p = row.r1_p, r1_c = row.r1_c, r1_w = hasw ? row.r1_w : nothing,
            p_res0_dual = row.pres0_d, p_res0_prim = row.pres0_p,
            c_res0_dual = row.cres0_d, c_res0_prim = row.cres0_p,
            w_res0_dual = hasw ? row.wres0_d : nothing, w_res0_prim = hasw ? row.wres0_p : nothing,
            p_res_exit = row.pres_exit, c_res_exit = row.cres_exit,
            w_res_exit = hasw ? row.wres_exit : nothing,
            Ldiag_min = minimum(abs, Ld), Ldiag_max = maximum(abs, Ld),
            ipm_pres = row.pres, ipm_dres = row.dres,
            bar_hdiag_med = row.bar_hdiag_med, bar_hdiag_frac_mid = row.bar_hdiag_frac_mid)
    # append the Gauss-quadrature harvest of the predictor rhs measure (nodes θ, weights w, len 10 NaN-padded)
    θ = row.ritz_θ_p; wq = row.ritz_w_p
    ritzcols = (; ritz_beta = row.ritz_beta_p, omm = row.omm_p,
        (Symbol(:ritz_t, j) => (j ≤ length(θ)  ? θ[j]  : T(NaN)) for j in 1:10)...,
        (Symbol(:ritz_w, j) => (j ≤ length(wq) ? wq[j] : T(NaN)) for j in 1:10)...)
    return merge(base, ritzcols, (chosen = false,))
end

# Serialize solve_logged records to a CSV at `path`. Generic over the record fields; enums are
# stringified, `nothing` becomes empty, floats use scientific notation (NaN preserved).
_csvcell(::Nothing) = ""
_csvcell(x::AbstractFloat) = isnan(x) ? "NaN" : @sprintf("%.15e", x)   # full precision: nodes near 1 carry info in last digits
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

function solve_logged(s0::AbstractSolver, grid = DEFAULT_ALPHA_GRID; itmax::Integer = s0.settings.itmax)
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
        for α in gs
            sc = deepcopy(s)
            sc.α[] = α
            st = step!(sc)
            push!(iter_recs, _oracle_record(i, α, sc, st))
            score = (iter_recs[end].state, iter_recs[end].ncraig)
            if score < bestscore                   # strict: first (lowest) α at the best score wins
                bestscore = score
                best = sc
                beststatus = st
                bestidx = lastindex(iter_recs)
            end
        end
        # skip the terminal-at-entry sweep: if the entry state is already converged every candidate
        # trivially returns ncraig=0 at a terminal status — record nothing and stop (removes the artifact).
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
