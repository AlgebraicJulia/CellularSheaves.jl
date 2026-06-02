const STALK_DIM = 2
const SMALL_SIZES = [20, 100]
const LARGE_SIZES = [500]
const ALL_SIZES = vcat(SMALL_SIZES, LARGE_SIZES)
const GRAPH_FAMILIES = ["cycle", "path"]

function make_path_sheaf(n, d)
    g = path_graph(n)
    sheaf_from_graph(g, d, k -> Matrix{Float64}(I, k, k); symmetric_edges=true)
end

function make_cycle_sheaf(n, d)
    g = cycle_graph(n)
    sheaf_from_graph(g, d, k -> Matrix{Float64}(I, k, k); symmetric_edges=true)
end

const FACTORIES = Dict(
    "cycle" => make_cycle_sheaf,
    "path" => make_path_sheaf,
)

function boundary_fixture(family, n)
    v2 = family == "path" ? n : n ÷ 2
    Dict(1 => zeros(STALK_DIM), v2 => ones(STALK_DIM))
end

function build_suite(sizes::Vector{Int}=ALL_SIZES)
    sheaves = Dict(
        (family, n) => FACTORIES[family](n, STALK_DIM)
        for family in GRAPH_FAMILIES
        for n in sizes
    )
    harmonic_boundaries = Dict(
        (family, n) => boundary_fixture(family, n)
        for family in GRAPH_FAMILIES
        for n in sizes
    )
    rng = MersenneTwister(42)
    global_section_inputs = Dict(
        (family, n) => randn(rng, STALK_DIM * n)
        for family in GRAPH_FAMILIES
        for n in sizes
    )

    suite = BenchmarkGroup()

    suite["coboundary_map"] = BenchmarkGroup()
    for family in GRAPH_FAMILIES
        suite["coboundary_map"][family] = BenchmarkGroup()
        for n in sizes
            s = sheaves[(family, n)]
            suite["coboundary_map"][family]["n$n"] = @benchmarkable coboundary_map($s)
        end
    end

    suite["laplacian"] = BenchmarkGroup()
    suite["laplacian"]["matrix_direct"] = BenchmarkGroup()
    for family in GRAPH_FAMILIES
        suite["laplacian"]["matrix_direct"][family] = BenchmarkGroup()
        for n in sizes
            s = sheaves[(family, n)]
            suite["laplacian"]["matrix_direct"][family]["n$n"] = @benchmarkable sheaf_laplacian_matrix_direct($s)
        end
    end

    suite["harmonic_extension"] = BenchmarkGroup()
    for family in GRAPH_FAMILIES
        suite["harmonic_extension"][family] = BenchmarkGroup()
        for n in sizes
            s = sheaves[(family, n)]
            boundary = harmonic_boundaries[(family, n)]
            suite["harmonic_extension"][family]["n$n"] = @benchmarkable harmonic_extension($s, $boundary)
        end
    end

    suite["nearest_global_section"] = BenchmarkGroup()
    suite["nearest_global_section"]["ldl"] = BenchmarkGroup()
    for family in GRAPH_FAMILIES
        suite["nearest_global_section"]["ldl"][family] = BenchmarkGroup()
        for n in sizes
            s = sheaves[(family, n)]
            x0 = global_section_inputs[(family, n)]
            suite["nearest_global_section"]["ldl"][family]["n$n"] = @benchmarkable nearest_global_section($s, $x0; method=:ldl)
        end
    end

    suite
end

benchmark_id(ids) = join(string.(ids), "/")

function suite_leaf_ids(group::BenchmarkGroup)
    sort!(String[benchmark_id(ids) for (ids, _) in BenchmarkTools.leaves(group)])
end

function _filter_group(group::BenchmarkGroup, prefix::Vector{String}, selected_ids::Set{String})
    filtered = BenchmarkGroup()
    child_count = 0
    for (key, value) in pairs(group)
        key_string = string(key)
        next_prefix = [prefix; key_string]
        if value isa BenchmarkGroup
            child = _filter_group(value, next_prefix, selected_ids)
            if child !== nothing
                filtered[key] = child
                child_count += 1
            end
        else
            leaf_id = benchmark_id(next_prefix)
            if leaf_id in selected_ids
                filtered[key] = value
                child_count += 1
            end
        end
    end
    child_count == 0 ? nothing : filtered
end

function filter_suite(group::BenchmarkGroup, selected_ids::Vector{String})
    filtered = _filter_group(group, String[], Set(selected_ids))
    filtered === nothing && error("The selected shard does not contain any benchmark leaves.")
    filtered
end
