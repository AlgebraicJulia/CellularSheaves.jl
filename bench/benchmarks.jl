using BenchmarkTools
using CellularSheaves
using Graphs
using LinearAlgebra
using Random

const SUITE = BenchmarkGroup()

# ---- Fixtures ----

function make_path_sheaf(n, d)
    g = path_graph(n)
    sheaf_from_graph(g, d, k -> Matrix{Float64}(I, k, k); symmetric_edges=true)
end

function make_cycle_sheaf(n, d)
    g = cycle_graph(n)
    sheaf_from_graph(g, d, k -> Matrix{Float64}(I, k, k); symmetric_edges=true)
end

const SIZES = [20, 100, 500]
const STALK_DIM = 2
const FACTORIES = [("path", make_path_sheaf), ("cycle", make_cycle_sheaf)]

# Pre-build all sheaves so fixture construction is not included in timings.
const SHEAVES = Dict(
    (name, n) => factory(n, STALK_DIM)
    for (name, factory) in FACTORIES
    for n in SIZES
)

# ---- coboundary_map ----

SUITE["coboundary_map"] = BenchmarkGroup()
for (name, _) in FACTORIES
    SUITE["coboundary_map"][name] = BenchmarkGroup()
    for n in SIZES
        s = SHEAVES[(name, n)]
        SUITE["coboundary_map"][name][n] = @benchmarkable coboundary_map($s)
    end
end

# ---- laplacian assembly ----

SUITE["laplacian"] = BenchmarkGroup()
for (name, _) in FACTORIES
    SUITE["laplacian"][name] = BenchmarkGroup()
    #SUITE["laplacian"][name]["matrix"] = BenchmarkGroup()
    SUITE["laplacian"][name]["matrix_direct"] = BenchmarkGroup()
    for n in SIZES
        s = SHEAVES[(name, n)]
        #SUITE["laplacian"][name]["matrix"][n] = @benchmarkable sheaf_laplacian_matrix($s)
        SUITE["laplacian"][name]["matrix_direct"][n] = @benchmarkable sheaf_laplacian_matrix_direct($s)
    end
end

# ---- harmonic extension ----
# path: fix both endpoints; cycle: fix vertex 1 and the midpoint vertex.

SUITE["harmonic_extension"] = BenchmarkGroup()
for (name, _) in FACTORIES
    SUITE["harmonic_extension"][name] = BenchmarkGroup()
    for n in SIZES
        s = SHEAVES[(name, n)]
        v2 = name == "path" ? n : n ÷ 2
        boundary = Dict(1 => zeros(STALK_DIM), v2 => ones(STALK_DIM))
        SUITE["harmonic_extension"][name][n] = @benchmarkable harmonic_extension($s, $boundary)
    end
end

# ---- nearest_global_section ----

SUITE["nearest_global_section"] = BenchmarkGroup()
Random.seed!(42)
for (name, _) in FACTORIES
    SUITE["nearest_global_section"][name] = BenchmarkGroup()
    for n in SIZES
        s = SHEAVES[(name, n)]
        x0 = randn(STALK_DIM * n)
        SUITE["nearest_global_section"][name][n] = @benchmarkable nearest_global_section($s, $x0; method=:ldl)
    end
end

# ---- Run ----

results = BenchmarkTools.run(SUITE, verbose=true; seconds=10, samples=100)

println("\n=== Median results ===")
display(median(results))

results_path = joinpath(@__DIR__, "results.json")
BenchmarkTools.save(results_path, results)
println("\nSaved $results_path")

# Compare against committed baseline, warn on any >2x regression.
baseline_path = joinpath(@__DIR__, "baseline.json")
if isfile(baseline_path)
    baseline = BenchmarkTools.load(baseline_path)[1]
    j = judge(median(results), median(baseline); time_tolerance=0.05)
    println("\n=== Regression report vs baseline ===")
    display(j)
    local any_regression = false
    for (keys, jmt) in leaves(j)
        if jmt.time === :regression && ratio(jmt).time > 2.0
            @warn "Regression >2x" benchmark = join(string.(keys), "/") factor = ratio(jmt).time
            any_regression = true
        end
    end
    any_regression || println("No regressions detected.")
else
    println("\nNo baseline.json found. To set one:")
    println("  cp bench/results.json bench/baseline.json")
end
