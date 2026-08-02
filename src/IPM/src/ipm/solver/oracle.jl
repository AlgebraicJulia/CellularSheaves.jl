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

# NOTE: the spectral harvest is no longer computed here. Every base solve (`solve_uzw!` in kkt/uzawa.jl,
# via `gk_spectral`) reads B's rhs measure AFTER craig!, on the exact rhs craig! consumed, at BOTH kmax=10
# and kmax=ncraig (a stand-in for a future inline-CRAIG harvest). Each of the step's base solves — predictor,
# corrector, and (HSD) woodbury — rides a SpectralHarvest per variant on the history row (ph10/phN/…/whN).

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

# Flatten one SpectralHarvest into CSV columns under a role/variant prefix (e.g. :p10, :pn, :c10, :cn,
# :w10, :wn): the σ̂² extremes, GK truncation β, 1-μ̂max, the actual GK step count k, and the top-10
# nodes/weights (NaN-padded).
function _harvest_cols(h::SpectralHarvest{T}, pre::Symbol) where {T}
    return (; Symbol(pre, :_s2min) => h.s2min,
              Symbol(pre, :_s2max) => h.s2max,
              Symbol(pre, :_beta)  => h.ritz_beta,
              Symbol(pre, :_omm)   => h.omm,
              Symbol(pre, :_k)     => h.k,
              (Symbol(pre, :_t, j) => h.θ[j] for j in 1:10)...,
              (Symbol(pre, :_w, j) => h.w[j] for j in 1:10)...)
end

# Flatten one PassTrace into CSV columns under a role prefix (:p, :c, :w): the per-pass ENTRY residual
# (dual/primal/τ split) and per-pass CRAIG count for refinement passes 1..PASSTRACE_LEN. NaN residual =
# the loop never reached that pass; nkry = −1 = that pass did not fire. See PassTrace docstring for the
# base-invocation (pass 0) mapping (r0_* / pbase−1 / this trace's pass-1 entry).
function _pass_cols(t::PassTrace{T}, pre::Symbol) where {T}
    return (; (Symbol(pre, :_dres, j) => t.dres[j] for j in 1:PASSTRACE_LEN)...,
              (Symbol(pre, :_pres, j) => t.pres[j] for j in 1:PASSTRACE_LEN)...,
              (Symbol(pre, :_tres, j) => t.tres[j] for j in 1:PASSTRACE_LEN)...,
              (Symbol(pre, :_nkry, j) => t.nkry[j] for j in 1:PASSTRACE_LEN)...)
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
            tau_stall = sc.settings.refine_stall, refine_cap = sc.settings.refine_itmax,   # refinement stall threshold + pass cap actually in force (run1 vs run2)
            normB2 = norm(sc.B)^2,                    # ordinary ‖B‖² (Frobenius) — the OTHER sense of "σ²max"; sigma2max below is the preconditioned Ritz σ̂²max
            sigma2min = row.ph10.s2min, sigma2max = row.ph10.s2max,   # predictor kmax=10 harvest (back-compat alias) — preconditioned B F⁻¹ Bᵀ Ritz, NOT ‖B‖²
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
    # append the two-variant (kmax=10 / kmax=ncraig) spectral harvests for every base solve this step —
    # predictor & corrector always, woodbury on HSD — each flattened under its role/variant prefix.
    ritzcols = merge(_harvest_cols(row.ph10, :p10), _harvest_cols(row.phN, :pn),
                     _harvest_cols(row.ch10, :c10), _harvest_cols(row.chN, :cn))
    if hasw
        ritzcols = merge(ritzcols, _harvest_cols(row.wh10, :w10), _harvest_cols(row.whN, :wn))
    end
    # per-pass refinement trace (entry residuals + CRAIG per pass) read from the workspace the captured
    # step just filled — predictor & corrector always, woodbury on HSD.
    passcols = merge(_pass_cols(sc.kkt.ptrace, :p), _pass_cols(sc.kkt.ctrace, :c))
    if hasw
        passcols = merge(passcols, _pass_cols(sc.kkt.wtrace, :w))
    end
    return merge(base, ritzcols, passcols, (chosen = false,))
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

############################################################################################
# solve_logged_fine — two-stage floor-resolving sweep (fine1 spec)
############################################################################################
#
# Per IPM iteration: (1) run the coarse half-decade grid exactly as solve_logged (same score, same
# advancement, so the iterate trajectory is identical); (2) locate the provisional floor from the coarse
# state==0 rows — nmin = min ncraig, α_floor = smallest α whose ncraig ≤ nmin+1; (3) re-run the frozen
# iterate's step on a `fstep`-decade lattice spanning [log α_floor ± halfwidth], skipping α within
# fstep/2 decades of a coarse point already run. Skip the fine pass (noting why) when there are no
# state==0 rows or the floor sits at the coarse-grid minimum (not bracketed). Every row carries a `grid`
# tag ("coarse"|"fine"). Returns (final_solver, records, notes) — notes is one NamedTuple per iteration
# for the summary (incl. a determinism self-check: α_floor re-solved once more, compared to its coarse row).
function solve_logged_fine(s0::AbstractSolver, coarse = DEFAULT_ALPHA_GRID;
                           itmax::Integer = s0.settings.itmax, halfwidth::Real = 1.0, fstep::Real = 0.1)
    s = s0
    records = NamedTuple[]
    notes = NamedTuple[]
    gc = sort(collect(coarse))
    logcoarse = log10.(gc)
    i = 0
    while i < itmax
        i += 1
        # ---- coarse pass (identical selection/advancement to solve_logged) ----
        iter_recs = NamedTuple[]
        best = nothing; bestscore = (typemax(Int), typemax(Int)); beststatus = CONTINUE; bestidx = 0
        for α in gc
            sc = deepcopy(s); sc.α[] = α; st = step!(sc)
            push!(iter_recs, merge(_oracle_record(i, α, sc, st), (grid = "coarse",)))
            score = (iter_recs[end].state, iter_recs[end].ncraig)
            if score < bestscore
                bestscore = score; best = sc; beststatus = st; bestidx = lastindex(iter_recs)
            end
        end
        if beststatus !== CONTINUE && iter_recs[bestidx].ncraig == 0
            break
        end
        iter_recs[bestidx] = merge(iter_recs[bestidx], (chosen = true,))
        # ---- locate provisional floor from coarse state==0 rows ----
        st0 = [r for r in iter_recs if r.state == 0]
        did_fine = false; reason = ""; αfl = NaN; nfine = 0; nmin = -1; detdiff = NaN
        if isempty(st0)
            reason = "no state==0 rows"
        else
            nmin = minimum(r.ncraig for r in st0)
            αfl = minimum(r.alpha for r in st0 if r.ncraig <= nmin + 1)
            if αfl <= first(gc) * (1 + 1e-9)
                reason = "floor at coarse-grid minimum (not bracketed)"
            else
                did_fine = true
                logfl = log10(αfl)
                for d in -halfwidth:fstep:halfwidth
                    lα = logfl + d
                    minimum(abs.(lα .- logcoarse)) < fstep / 2 && continue   # skip coincident-with-coarse
                    α = exp10(lα)
                    sc = deepcopy(s); sc.α[] = α; st = step!(sc)
                    push!(iter_recs, merge(_oracle_record(i, α, sc, st), (grid = "fine",)))
                    nfine += 1
                end
                # determinism self-check: re-solve α_floor once more on the same frozen iterate.
                crow = iter_recs[findfirst(r -> r.alpha == αfl, iter_recs)]
                sc = deepcopy(s); sc.α[] = αfl; st = step!(sc)
                rr = _oracle_record(i, αfl, sc, st)
                rel(a, b) = abs(a - b) / max(abs(b), one(a))
                detdiff = max(rel(rr.ncraig, crow.ncraig), rel(rr.pbase, crow.pbase), rel(rr.r0_p, crow.r0_p))
            end
        end
        push!(notes, (iter = i, alpha_floor = αfl, nmin = nmin, did_fine = did_fine,
                      nfine = nfine, det_reldiff = detdiff, reason = reason))
        append!(records, iter_recs)
        s = best
        beststatus === CONTINUE || break
    end
    return s, records, notes
end

# 3-tread bracketing test: ≥2 fine rows at plateau (min ncraig), nmin+1, and nmin+2.
function _brackets3(frows)
    isempty(frows) && return false
    plat = minimum(r.ncraig for r in frows)
    return count(r -> r.ncraig == plat,     frows) >= 2 &&
           count(r -> r.ncraig == plat + 1, frows) >= 2 &&
           count(r -> r.ncraig >= plat + 2, frows) >= 2
end

# one edge's fine lattice at a given half-width (0.1-decade), coincident-with-coarse points skipped.
function _fine_lattice(s, i, edgeα, edgename, width, fstep, logcoarse)
    rows = NamedTuple[]
    logc = log10(edgeα)
    for d in -width:fstep:width
        lα = logc + d
        minimum(abs.(lα .- logcoarse)) < fstep / 2 && continue
        α = exp10(lα)
        sc = deepcopy(s); sc.α[] = α; st = step!(sc)
        push!(rows, merge(_oracle_record(i, α, sc, st), (grid = "fine", edge = edgename)))
    end
    return rows
end

# solve_logged_fine2 (fine2 spec) — like solve_logged_fine but brackets BOTH window edges (floor and
# ceiling, the two ends of the ncraig plateau: smallest/largest α with ncraig ≤ nmin+1) and widens each
# edge's lattice adaptively (widths ±1→1.5→2.5→3.5) until the 3-tread bracketing test passes, rather than
# skipping. Every row carries `grid` (coarse|fine) and `edge` (coarse|floor|ceil).
function solve_logged_fine2(s0::AbstractSolver, coarse = DEFAULT_ALPHA_GRID;
                            itmax::Integer = s0.settings.itmax,
                            widths = (1.0, 1.5, 2.5, 3.5), fstep::Real = 0.1)
    s = s0
    records = NamedTuple[]
    notes = NamedTuple[]
    gc = sort(collect(coarse)); logcoarse = log10.(gc)
    i = 0
    while i < itmax
        i += 1
        iter_recs = NamedTuple[]
        best = nothing; bestscore = (typemax(Int), typemax(Int)); beststatus = CONTINUE; bestidx = 0
        for α in gc
            sc = deepcopy(s); sc.α[] = α; st = step!(sc)
            push!(iter_recs, merge(_oracle_record(i, α, sc, st), (grid = "coarse", edge = "coarse")))
            score = (iter_recs[end].state, iter_recs[end].ncraig)
            if score < bestscore
                bestscore = score; best = sc; beststatus = st; bestidx = lastindex(iter_recs)
            end
        end
        if beststatus !== CONTINUE && iter_recs[bestidx].ncraig == 0
            break
        end
        iter_recs[bestidx] = merge(iter_recs[bestidx], (chosen = true,))
        st0 = [r for r in iter_recs if r.state == 0]
        if isempty(st0)
            push!(notes, (iter = i, edge = "both", did_fine = false, width = 0.0,
                          bracketed = false, nfine = 0, det_reldiff = NaN, reason = "no state==0 rows"))
        else
            nmin = minimum(r.ncraig for r in st0)
            cand = [r for r in st0 if r.ncraig <= nmin + 1]
            for (edgename, edgeα) in (("floor", minimum(r.alpha for r in cand)),
                                      ("ceil",  maximum(r.alpha for r in cand)))
                fine = NamedTuple[]; usedw = first(widths); brk = false
                for w in widths
                    fine = _fine_lattice(s, i, edgeα, edgename, w, fstep, logcoarse)
                    usedw = w
                    _brackets3(fine) && (brk = true; break)
                end
                crow = iter_recs[findfirst(r -> r.alpha == edgeα, iter_recs)]
                sc = deepcopy(s); sc.α[] = edgeα; st = step!(sc)
                rr = _oracle_record(i, edgeα, sc, st)
                rel(a, b) = abs(a - b) / max(abs(b), one(a))
                dd = max(rel(rr.ncraig, crow.ncraig), rel(rr.pbase, crow.pbase), rel(rr.r0_p, crow.r0_p))
                append!(iter_recs, fine)
                push!(notes, (iter = i, edge = edgename, did_fine = true, width = usedw,
                              bracketed = brk, nfine = length(fine), det_reldiff = dd,
                              reason = brk ? "" : "not bracketed at max width"))
            end
        end
        append!(records, iter_recs)
        s = best
        beststatus === CONTINUE || break
    end
    return s, records, notes
end

############################################################################################
# snapshot_capture / oracle_mu_trajectory — pre-solve KKT input state capture (snapshots spec)
############################################################################################
#
# Advance along the oracle's chosen-α trajectory (same selection as solve_logged) to the ENTRY of
# iteration K, then run the chosen-α step once more with the SNAP hook armed, recording every base
# solve's pre-solve input state (A, B, f, g, y0). No solver logic changes: capture is a guarded
# observer, off in production. Returns nothing if the solve terminates before reaching K.
function snapshot_capture(s0::AbstractSolver, K::Integer, grid = DEFAULT_ALPHA_GRID)
    s = s0
    gs = sort(collect(grid))
    chosenα = NaN
    i = 0
    while i < K
        i += 1
        best = nothing; bestscore = (typemax(Int), typemax(Int)); bestα = NaN; beststatus = CONTINUE
        for α in gs
            sc = deepcopy(s); sc.α[] = α; st = step!(sc); row = sc.hist[end]
            score = (_oracle_state(row), _oracle_craig(row))
            if score < bestscore
                bestscore = score; best = sc; bestα = α; beststatus = st
            end
        end
        best === nothing && return nothing
        if i == K
            chosenα = bestα
            break
        end
        s = best
        beststatus === CONTINUE || return nothing
    end
    # entry-state scalars (before the captured step), then the armed capture step at the chosen α
    entry = (mu = mu(s), nc = s.nc[], ng = s.ng[],
             tau = hasproperty(s, :τ) ? s.τ[] : NaN, kappa = hasproperty(s, :κ) ? s.κ[] : NaN)
    SNAP[] = NamedTuple[]
    sc = deepcopy(s); sc.α[] = chosenα; st = step!(sc)
    snaps = SNAP[]::Vector; SNAP[] = nothing
    row = sc.hist[end]
    T = eltype(sc.wrk.rp)
    μ1 = first(sc.hist.μ)
    ftol  = min(sc.settings.forcing_frac * row.μ / μ1, sc.settings.forcing_ceil)
    fltol = 100 * eps(T) * (one(T) + max(norm(sc.wrk.rp, Inf), norm(sc.wrk.rd, Inf)))
    return (; chosenα, entry, snaps, row, status = st, force_tol = ftol, floor_tol = fltol)
end

# snapshot_capture_at (snapshots3 spec) — like snapshot_capture, but run the captured step at a SPECIFIED
# α (not the chosen one), so the same iterate can be captured on either side of the ρ-shift transition.
# The trajectory to iteration K's entry still follows the chosen-α path. Also reports norm_dp/norm_dy and
# the step's ρ (row.ρ), for the shifted-regime meta fields.
function snapshot_capture_at(s0::AbstractSolver, K::Integer, alpha::Real, grid = DEFAULT_ALPHA_GRID)
    s = s0
    gs = sort(collect(grid))
    i = 0
    while i < K
        i += 1
        best = nothing; bestscore = (typemax(Int), typemax(Int)); beststatus = CONTINUE
        for α in gs
            sc = deepcopy(s); sc.α[] = α; st = step!(sc); row = sc.hist[end]
            score = (_oracle_state(row), _oracle_craig(row))
            if score < bestscore
                bestscore = score; best = sc; beststatus = st
            end
        end
        best === nothing && return nothing
        i == K && break                       # reached iteration K's entry iterate (s)
        s = best
        beststatus === CONTINUE || return nothing
    end
    entry = (mu = mu(s), nc = s.nc[], ng = s.ng[],
             tau = hasproperty(s, :τ) ? s.τ[] : NaN, kappa = hasproperty(s, :κ) ? s.κ[] : NaN)
    SNAP[] = NamedTuple[]
    sc = deepcopy(s); sc.α[] = oftype(sc.α[], alpha); st = step!(sc)
    snaps = SNAP[]::Vector; SNAP[] = nothing
    row = sc.hist[end]
    T = eltype(sc.wrk.rp)
    μ1 = first(sc.hist.μ)
    ftol  = min(sc.settings.forcing_frac * row.μ / μ1, sc.settings.forcing_ceil)
    fltol = 100 * eps(T) * (one(T) + max(norm(sc.wrk.rp, Inf), norm(sc.wrk.rd, Inf)))
    return (; alpha = sc.α[], entry, snaps, row, status = st, force_tol = ftol, floor_tol = fltol,
            norm_dp = norm(sc.wrk.Δp), norm_dy = norm(sc.wrk.Δy))
end

# μ per iteration along the chosen-α trajectory (for matching the X03 ipm/hsd pair by μ).
function oracle_mu_trajectory(s0::AbstractSolver, maxit::Integer, grid = DEFAULT_ALPHA_GRID)
    s = s0; gs = sort(collect(grid)); mus = Float64[]
    for _ in 1:maxit
        best = nothing; bestscore = (typemax(Int), typemax(Int)); beststatus = CONTINUE
        for α in gs
            sc = deepcopy(s); sc.α[] = α; st = step!(sc); row = sc.hist[end]
            score = (_oracle_state(row), _oracle_craig(row))
            if score < bestscore
                bestscore = score; best = sc; beststatus = st
            end
        end
        best === nothing && break
        push!(mus, best.hist[end].μ)
        s = best
        beststatus === CONTINUE || break
    end
    return mus
end
