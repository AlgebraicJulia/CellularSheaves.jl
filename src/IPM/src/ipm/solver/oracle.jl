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

const DEFAULT_ALPHA_GRID = [round(10.0^e, sigdigits=4) for e in 0.0:0.5:10.0]   # half-decades 1e0..1e10

# σ̂²_min from a Golub-Kahan lower-bidiagonal (dv=diagonal αₖ, ev=subdiagonal βₖ): μ = smallest Ritz
# value = min singular value², recovered to the augmented σ² via μ/(α(1-μ)) (μ∈(0,1) for F=βH+BᵀB).
function ritz_sigma2min(dv::AbstractVector{T}, ev::AbstractVector{T}, α_aug::T) where {T}
    k = length(dv)
    k ≥ 1 || return T(NaN)
    Bk = k == 1 ? Bidiagonal(T[dv[1]], T[], :L) : Bidiagonal(collect(dv), collect(@view ev[1:k-1]), :L)
    μ = minimum(svdvals(Bk))^2
    return (zero(T) < μ < one(T)) ? μ / (α_aug * (one(T) - μ)) : T(NaN)
end

# Standalone σ̂²_min estimate — DECOUPLED from the base solve. Runs its own k-step preconditioned
# Golub-Kahan bidiagonalization of B under N = F⁻¹ (reusing only the factorization F/divwrk), with
# internally-allocated work vectors. Nothing in the live solve is read or mutated; the extra ldivs
# are off-books (the oracle scores the real step!). Returns σ̂²_min (NaN if degenerate).
function oracle_sigma2min(kkt, B, α_aug::T; kmax::Int = 10, b::Union{Nothing, AbstractVector{T}} = nothing) where {T}
    m, n = size(B)
    F = kkt.F; dw = kkt.divwrk
    u = b === nothing ? ones(T, m) : copy(b)
    Nv = zeros(T, n); v = zeros(T, n); Atu = zeros(T, n); Av = zeros(T, m)
    dv = T[]; ev = T[]
    β = norm(u); β == zero(T) && return T(NaN)
    u ./= β
    for _ in 1:kmax
        mul!(Atu, B', u)                                 # Bᵀ u
        @. Nv = Atu - β * Nv
        copyto!(v, Nv); ldiv!(dw, F, v); ldiv!(dw, F', v)   # v = F⁻¹ Nv
        α = sqrt(max(real(dot(v, Nv)), zero(T)))            # elliptic norm √(v·Nv)
        α == zero(T) && break
        push!(dv, α)
        v ./= α; Nv ./= α
        mul!(Av, B, v)                                    # B v
        @. u = Av - α * u
        β = norm(u)
        push!(ev, β)
        β == zero(T) && break
        u ./= β
    end
    return ritz_sigma2min(dv, ev, α_aug)
end

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
    return (iter = i, alpha = α, state = _oracle_state(row), ncraig = _oracle_craig(row),
            ipm_status = st, mu = μ, mu_next = mu(sc), rho = row.ρ,
            force_tol = ftol, floor_tol = fltol, sigma2min = oracle_sigma2min(sc.kkt, sc.B, T(α)),
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
            chosen = false)
end

# Serialize solve_logged records to a CSV at `path`. Generic over the record fields; enums are
# stringified, `nothing` becomes empty, floats use scientific notation (NaN preserved).
_csvcell(::Nothing) = ""
_csvcell(x::AbstractFloat) = isnan(x) ? "NaN" : @sprintf("%.6e", x)
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
