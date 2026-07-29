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

function _oracle_record(i, α, row, st)
    return (iter = i, alpha = α, state = _oracle_state(row), ncraig = _oracle_craig(row),
            ipm_status = st, pstat = row.pstat, cstat = row.cstat,
            wstat = hasproperty(row, :wstat) ? row.wstat : nothing,
            pbase = row.pbase, prefn = row.prefn, cbase = row.cbase, crefn = row.crefn,
            wbase = hasproperty(row, :wbase) ? row.wbase : 0,
            wrefn = hasproperty(row, :wrefn) ? row.wrefn : 0,
            chosen = false)
end

function solve_logged(s0::AbstractSolver, grid = DEFAULT_ALPHA_GRID; itmax::Integer = s0.settings.itmax)
    s = s0
    records = NamedTuple[]
    gs = sort(collect(grid))                       # ascending ⇒ ties on (state,ncraig) → lowest α
    i = 0
    while i < itmax
        i += 1
        best = nothing
        bestscore = (typemax(Int), typemax(Int))
        beststatus = CONTINUE
        bestidx = 0
        for α in gs
            sc = deepcopy(s)
            sc.α[] = α
            st = step!(sc)
            row = sc.hist[end]
            push!(records, _oracle_record(i, α, row, st))
            score = (_oracle_state(row), _oracle_craig(row))
            if score < bestscore                   # strict: first (lowest) α at the best score wins
                bestscore = score
                best = sc
                beststatus = st
                bestidx = lastindex(records)
            end
        end
        records[bestidx] = merge(records[bestidx], (chosen = true,))
        s = best
        beststatus === CONTINUE || break           # terminal status at the chosen α ⇒ done
    end
    return s, records
end
