const SHARD_ORDER = [
    "assembly-small",
    "solver-small",
    "extension-small",
    "assembly-large",
    "solver-large",
    "extension-large",
]

function benchmark_ids(operation::String, sizes::Vector{Int})
    if operation == "coboundary_map"
        return [join((operation, family, "n$n"), "/") for family in GRAPH_FAMILIES for n in sizes]
    elseif operation == "laplacian"
        return [join((operation, "matrix_direct", family, "n$n"), "/") for family in GRAPH_FAMILIES for n in sizes]
    elseif operation == "harmonic_extension"
        return [join((operation, family, "n$n"), "/") for family in GRAPH_FAMILIES for n in sizes]
    elseif operation == "nearest_global_section"
        return [join((operation, "ldl", family, "n$n"), "/") for family in GRAPH_FAMILIES for n in sizes]
    end
    error("Unknown benchmark operation: $operation")
end

const SHARD_MANIFEST = Dict(
    "assembly-small" => vcat(benchmark_ids("coboundary_map", SMALL_SIZES), benchmark_ids("laplacian", SMALL_SIZES)),
    "solver-small" => benchmark_ids("nearest_global_section", SMALL_SIZES),
    "extension-small" => benchmark_ids("harmonic_extension", SMALL_SIZES),
    "assembly-large" => vcat(benchmark_ids("coboundary_map", LARGE_SIZES), benchmark_ids("laplacian", LARGE_SIZES)),
    "solver-large" => benchmark_ids("nearest_global_section", LARGE_SIZES),
    "extension-large" => benchmark_ids("harmonic_extension", LARGE_SIZES),
)

available_profiles() = ["small", "large", "full"]

function normalize_profile(profile::AbstractString)
    profile in available_profiles() || error("Unknown benchmark profile '$profile'. Expected one of $(join(available_profiles(), ", ")).")
    profile
end

function available_shards(profile::AbstractString)
    normalized = normalize_profile(profile)
    if normalized == "small"
        return filter(shard -> endswith(shard, "-small"), SHARD_ORDER)
    elseif normalized == "large"
        return filter(shard -> endswith(shard, "-large"), SHARD_ORDER)
    end
    copy(SHARD_ORDER)
end

function runner_kind_from_env()
    if haskey(ENV, "GITHUB_ACTIONS")
        return "github-actions"
    elseif haskey(ENV, "SLURM_JOB_ID")
        return "slurm"
    end
    "local"
end

function validate_profile_for_runner(profile::AbstractString)
    runner_kind = runner_kind_from_env()
    if runner_kind == "github-actions" && profile != "small" && get(ENV, "BENCHMARK_ALLOW_LARGE_IN_CI", "false") != "true"
        error("GitHub Actions runs are restricted to BENCHMARK_PROFILE=small unless BENCHMARK_ALLOW_LARGE_IN_CI=true is set.")
    end
    nothing
end

function selected_shards(profile::AbstractString, shard::AbstractString)
    normalized = normalize_profile(profile)
    validate_profile_for_runner(normalized)
    if shard == "all"
        return available_shards(normalized)
    end
    shard in keys(SHARD_MANIFEST) || error("Unknown benchmark shard '$shard'.")
    shard in available_shards(normalized) || error("Shard '$shard' is not part of BENCHMARK_PROFILE=$normalized.")
    [shard]
end

function selected_leaf_ids(profile::AbstractString, shard::AbstractString)
    shard_names = selected_shards(profile, shard)
    sort!(unique!(vcat([copy(SHARD_MANIFEST[name]) for name in shard_names]...)))
end

function filtered_suite(profile::AbstractString, shard::AbstractString)
    filter_suite(build_suite(), selected_leaf_ids(profile, shard))
end

function filtered_suite_from_env()
    filtered_suite(get(ENV, "BENCHMARK_PROFILE", "small"), get(ENV, "BENCHMARK_SHARD", "all"))
end
