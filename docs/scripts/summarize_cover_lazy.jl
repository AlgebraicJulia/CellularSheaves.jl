# Summary of tropical_cover_lazy.jl runs on one map, under the MAPF benchmark protocol
# (per scenario, the first N agents for growing N; Stern et al. 2019), at one or more
# time limits. Reads docs/figures/tropical_cover_lazy_<map>_s<k>-<k>.csv.
#
#   julia docs/scripts/summarize_cover_lazy.jl MAP [limits...]     (default limits 30 300)
# SUMMARY_BINS=1,11,21,31,51 sets the first N of each bin (default: bins of 10, or of 100).
#
# Per scenario, at each time limit (cumulative time to reach N, i.e. from scratch):
#   - planned: the largest N with a prioritized plan (the upper bound);
#   - certified: the largest N up to which every instance is certified optimal (the
#     cover bound meets the plan);
# and, over all scenarios and instances within the limit, binned by N:
#   - the share of instances certified;
#   - the share of the gap between the decoupled bound and the plan closed by the
#     cover of edges, and by the cover with triangles;
#   - triangles evaluated and obstructed (o(σ) > 0).

using Printf
using Statistics

mapname = ARGS[1]
limits = length(ARGS) >= 2 ? parse.(Float64, ARGS[2:end]) : [30.0, 300.0]
const BIN_STARTS = parse.(Int, filter(!isempty, split(get(ENV, "SUMMARY_BINS", ""), ",")))   # e.g. 1,11,21,31,51
dir = joinpath(@__DIR__, "..", "figures")
const TAG = get(ENV, "SUMMARY_TAG", "")                            # e.g. "_opt" for the runs with optima
files = filter(f -> occursin(Regex("^tropical_cover_lazy_\\Q$(mapname)\\E_s\\d+-\\d+\\Q$(TAG)\\E\\.csv\$"), f), readdir(dir))
rows = NamedTuple[]
for f in files
    lines = readlines(joinpath(dir, f))
    header = Symbol.(split(lines[1], ","))
    for l in lines[2:end]
        v = split(l, ",")
        parsev(s) = s == "true" ? true : s == "false" ? false : something(tryparse(Int, s), tryparse(Float64, s), s)
        push!(rows, NamedTuple{Tuple(header)}(Tuple(parsev.(v))))
    end
end
scenarios = sort(unique(r.scenario for r in rows))
@printf("%s: %d scenarios, %d instances\n", mapname, length(scenarios), length(rows))

out = joinpath(dir, "tropical_cover_lazy_$(mapname)$(TAG)_summary.txt")
open(out, "w") do io
    for io2 in (stdout, io)
        @printf(io2, "%s: %d scenario files, %d instances\n", mapname, length(scenarios), length(rows))
    end
    for T in limits
        planned, certified = Int[], Int[]
        for k in scenarios
            rs = sort(filter(r -> r.scenario == k && r.seconds <= T, rows); by = r -> r.N)
            push!(planned, maximum((r.N for r in rs if isfinite(r.upper)); init = 0))
            c = 0
            for r in rs
                r.certified || break
                c = r.N
            end
            push!(certified, c)
        end
        within = filter(r -> r.seconds <= T && isfinite(r.upper), rows)
        lines = String[]
        push!(lines, @sprintf("\ntime limit %.0f s", T))
        push!(lines, @sprintf("  largest N with a plan, per scenario: mean %.1f, median %.0f, range %d–%d",
            mean(planned), median(planned), minimum(planned), maximum(planned)))
        push!(lines, @sprintf("  largest N certified optimal throughout: mean %.1f, median %.0f, range %d–%d",
            mean(certified), median(certified), minimum(certified), maximum(certified)))
        Nmax = maximum(r.N for r in within)
        width = Nmax <= 100 ? 10 : 100
        starts = isempty(BIN_STARTS) ? collect(1:width:Nmax) : BIN_STARTS
        push!(lines, "  N bin      instances  certified  gap closed: edges / +triangles   triangles: evaluated / obstructed   packed")
        for (b, lo) in enumerate(starts)
            hi = b < length(starts) ? starts[b + 1] - 1 : max(Nmax, lo)
            rs = filter(r -> lo <= r.N <= hi, within)
            isempty(rs) && continue
            gap = sum(r.upper - r.decoupled for r in rs)
            cov(f) = gap == 0 ? "   –" : @sprintf("%3.0f%%", 100 * sum(f(r) - r.decoupled for r in rs) / gap)
            push!(lines, @sprintf("  %4d–%-4d  %9d  %8.0f%%  %13s / %s   %20.0f / %.1f   %5.0f%%", lo, hi, length(rs),
                100 * count(r -> r.certified, rs) / length(rs), cov(r -> r.bound_edges), cov(r -> r.bound_sheaf),
                mean(r.triangles_evaluated for r in rs), mean(r.triangles_obstructed for r in rs),
                100 * count(r -> r.packed_components > 0, rs) / length(rs)))
        end
        # Where the optimum is known (LAZY_OPT runs): the share of the true coordination
        # cost, optimum − decoupled, that the covers recover, and how often the plan is optimal.
        known = filter(r -> r.optimal isa Real && isfinite(r.optimal), within)
        if !isempty(known)
            push!(lines, "  with the optimum known:  instances  plan optimal  coordination cost recovered: edges / +triangles")
            for (b, lo) in enumerate(starts)
                hi = b < length(starts) ? starts[b + 1] - 1 : max(Nmax, lo)
                rs = filter(r -> lo <= r.N <= hi, known)
                isempty(rs) && continue
                cost = sum(r.optimal - r.decoupled for r in rs)
                rec(f) = cost == 0 ? "   –" : @sprintf("%3.0f%%", 100 * sum(f(r) - r.decoupled for r in rs) / cost)
                push!(lines, @sprintf("  %4d–%-4d  %20d  %11.0f%%  %24s / %s", lo, hi, length(rs),
                    100 * count(r -> r.upper == r.optimal, rs) / length(rs), rec(r -> r.bound_edges), rec(r -> r.bound_sheaf)))
            end
        end
        push!(lines, @sprintf("  capped pair searches: %d; capped triangle searches: %d; skipped triangles: %d; packed components: %d (instances)",
            sum(r.pairs_capped for r in within), sum(r.triangles_capped for r in within),
            sum(r.triangles_skipped for r in within), count(r -> r.packed_components > 0, within)))
        for io2 in (stdout, io), l in lines
            println(io2, l)
        end
    end
end
