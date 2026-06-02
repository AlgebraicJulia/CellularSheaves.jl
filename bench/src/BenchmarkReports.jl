escape_html(value) = replace(string(value), '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;", '\'' => "&#39;")

function format_seconds_ns(value)
    if value >= 1.0e9
        return @sprintf("%.3f s", value / 1.0e9)
    elseif value >= 1.0e6
        return @sprintf("%.3f ms", value / 1.0e6)
    elseif value >= 1.0e3
        return @sprintf("%.3f us", value / 1.0e3)
    end
    @sprintf("%.0f ns", value)
end

function format_bytes(value)
    if value >= 1024^2
        return @sprintf("%.2f MiB", value / 1024^2)
    elseif value >= 1024
        return @sprintf("%.2f KiB", value / 1024)
    end
    string(value, " B")
end

format_ratio(value) = value === nothing ? "n/a" : @sprintf("%.3fx", float(value))

function parse_bool_env(name::AbstractString, default::Bool=false)
    raw = lowercase(get(ENV, name, default ? "true" : "false"))
    raw in ("1", "true", "yes", "on")
end

function parse_int_env(name::AbstractString, default::Int)
    parse(Int, get(ENV, name, string(default)))
end

function parse_float_env(name::AbstractString, default::Float64)
    parse(Float64, get(ENV, name, string(default)))
end

function default_result_dir()
    joinpath(@__DIR__, "..", "results")
end

function ensure_dir(path::AbstractString)
    mkpath(path)
    path
end

function repo_root()
    normpath(joinpath(@__DIR__, "..", ".."))
end

function benchmark_metadata(profile::AbstractString, shard::AbstractString)
    Dict(
        "profile" => profile,
        "shard" => shard,
        "git_ref" => get(ENV, "BENCHMARK_REF", "working-tree"),
        "runner_kind" => runner_kind_from_env(),
        "hostname" => gethostname(),
        "timestamp" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
    )
end

function write_json(path::AbstractString, data)
    open(path, "w") do io
        JSON.print(io, data, 2)
    end
    path
end

function trial_row(benchmark_id::AbstractString, estimate, metadata::Dict{String,Any})
    Dict(
        "benchmark_id" => benchmark_id,
        "shard" => metadata["shard"],
        "profile" => metadata["profile"],
        "time_ns_median" => Int(estimate.time),
        "gc_time_ns_median" => Int(estimate.gctime),
        "memory_bytes" => Int(estimate.memory),
        "allocs" => Int(estimate.allocs),
        "git_ref" => metadata["git_ref"],
        "runner_kind" => metadata["runner_kind"],
        "hostname" => metadata["hostname"],
        "julia_version" => metadata["julia_version"],
        "threads" => metadata["threads"],
        "timestamp" => metadata["timestamp"],
    )
end

function summary_rows(group::BenchmarkGroup, metadata::Dict{String,Any})
    median_group = median(group)
    rows = Dict{String,Any}[]
    for (ids, estimate) in BenchmarkTools.leaves(median_group)
        push!(rows, trial_row(benchmark_id(ids), estimate, metadata))
    end
    sort!(rows; by=row -> row["benchmark_id"])
end

function summary_markdown(rows::Vector{Dict{String,Any}}, metadata::Dict{String,Any})
    io = IOBuffer()
    println(io, "# Benchmark shard report")
    println(io)
    println(io, "- **Shard:** `", metadata["shard"], "`")
    println(io, "- **Profile:** `", metadata["profile"], "`")
    println(io, "- **Runner:** `", metadata["runner_kind"], "` on `", metadata["hostname"], "`")
    println(io, "- **Git ref:** `", metadata["git_ref"], "`")
    println(io, "- **Julia:** `", metadata["julia_version"], "` with `", metadata["threads"], "` threads")
    println(io)
    println(io, "| Benchmark ID | Median time | Memory | Allocs |")
    println(io, "| --- | ---: | ---: | ---: |")
    for row in rows
        println(
            io,
            "| `", row["benchmark_id"], "` | ",
            format_seconds_ns(row["time_ns_median"]), " | ",
            format_bytes(row["memory_bytes"]), " | ",
            row["allocs"], " |",
        )
    end
    String(take!(io))
end

function comparison_row(benchmark_id::AbstractString, target_estimate, baseline_estimate, judgement, metadata::Dict{String,Any})
    time_ratio = baseline_estimate.time == 0 ? nothing : target_estimate.time / baseline_estimate.time
    memory_ratio = baseline_estimate.memory == 0 ? nothing : target_estimate.memory / baseline_estimate.memory
    label = BenchmarkTools.isregression(judgement) ? "regression" : BenchmarkTools.isimprovement(judgement) ? "improvement" : "invariant"
    Dict(
        "benchmark_id" => benchmark_id,
        "shard" => metadata["shard"],
        "profile" => metadata["profile"],
        "target_time_ns_median" => Int(target_estimate.time),
        "baseline_time_ns_median" => Int(baseline_estimate.time),
        "time_ratio" => time_ratio,
        "target_memory_bytes" => Int(target_estimate.memory),
        "baseline_memory_bytes" => Int(baseline_estimate.memory),
        "memory_ratio" => memory_ratio,
        "judgement" => label,
        "target_ref" => metadata["target_ref"],
        "baseline_ref" => metadata["baseline_ref"],
        "runner_kind" => metadata["runner_kind"],
        "hostname" => metadata["hostname"],
        "timestamp" => metadata["timestamp"],
    )
end

function comparison_rows(target_group::BenchmarkGroup, baseline_group::BenchmarkGroup, judgement_group::BenchmarkGroup, metadata::Dict{String,Any})
    target_map = Dict(benchmark_id(ids) => estimate for (ids, estimate) in BenchmarkTools.leaves(median(target_group)))
    baseline_map = Dict(benchmark_id(ids) => estimate for (ids, estimate) in BenchmarkTools.leaves(median(baseline_group)))
    judgement_map = Dict(benchmark_id(ids) => estimate for (ids, estimate) in BenchmarkTools.leaves(judgement_group))
    rows = Dict{String,Any}[]
    for id in sort(collect(keys(judgement_map)))
        push!(rows, comparison_row(id, target_map[id], baseline_map[id], judgement_map[id], metadata))
    end
    rows
end

function comparison_markdown(rows::Vector{Dict{String,Any}}, metadata::Dict{String,Any})
    io = IOBuffer()
    println(io, "# Benchmark comparison report")
    println(io)
    println(io, "- **Shard:** `", metadata["shard"], "`")
    println(io, "- **Profile:** `", metadata["profile"], "`")
    println(io, "- **Target ref:** `", metadata["target_ref"], "`")
    println(io, "- **Baseline ref:** `", metadata["baseline_ref"], "`")
    println(io)
    println(io, "| Benchmark ID | Target time | Baseline time | Time ratio | Target memory | Baseline memory | Memory ratio | Judgement |")
    println(io, "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    for row in rows
        println(
            io,
            "| `", row["benchmark_id"], "` | ",
            format_seconds_ns(row["target_time_ns_median"]), " | ",
            format_seconds_ns(row["baseline_time_ns_median"]), " | ",
            format_ratio(row["time_ratio"]), " | ",
            format_bytes(row["target_memory_bytes"]), " | ",
            format_bytes(row["baseline_memory_bytes"]), " | ",
            format_ratio(row["memory_ratio"]), " | `", row["judgement"], "` |",
        )
    end
    String(take!(io))
end

function combined_summary_markdown(rows::Vector{Dict{String,Any}}, expected_shards::Vector{String}, completed_shards::Vector{String}, profile::AbstractString)
    missing = sort(setdiff(expected_shards, completed_shards))
    io = IOBuffer()
    println(io, "# Benchmark report")
    println(io)
    println(io, "- **Profile:** `", profile, "`")
    println(io, "- **Expected shards:** `", join(expected_shards, ", "), "`")
    println(io, "- **Completed shards:** `", join(completed_shards, ", "), "`")
    println(io, "- **Missing shards:** `", isempty(missing) ? "none" : join(missing, ", "), "`")
    println(io)
    println(io, "| Benchmark ID | Shard | Median time | Memory | Allocs | Runner | Host |")
    println(io, "| --- | --- | ---: | ---: | ---: | --- | --- |")
    for row in rows
        println(
            io,
            "| `", row["benchmark_id"], "` | `", row["shard"], "` | ",
            format_seconds_ns(row["time_ns_median"]), " | ",
            format_bytes(row["memory_bytes"]), " | ",
            row["allocs"], " | `", row["runner_kind"], "` | `", row["hostname"], "` |",
        )
    end
    String(take!(io))
end

function combined_comparison_markdown(rows::Vector{Dict{String,Any}}, expected_shards::Vector{String}, completed_shards::Vector{String}, profile::AbstractString)
    missing = sort(setdiff(expected_shards, completed_shards))
    io = IOBuffer()
    println(io, "# Benchmark comparison report")
    println(io)
    println(io, "- **Profile:** `", profile, "`")
    println(io, "- **Expected shards:** `", join(expected_shards, ", "), "`")
    println(io, "- **Completed shards:** `", join(completed_shards, ", "), "`")
    println(io, "- **Missing shards:** `", isempty(missing) ? "none" : join(missing, ", "), "`")
    println(io)
    println(io, "| Benchmark ID | Shard | Target time | Baseline time | Time ratio | Judgement |")
    println(io, "| --- | --- | ---: | ---: | ---: | --- |")
    for row in rows
        println(
            io,
            "| `", row["benchmark_id"], "` | `", row["shard"], "` | ",
            format_seconds_ns(row["target_time_ns_median"]), " | ",
            format_seconds_ns(row["baseline_time_ns_median"]), " | ",
            format_ratio(row["time_ratio"]), " | `", row["judgement"], "` |",
        )
    end
    String(take!(io))
end

function report_html(title::AbstractString, table_headers::Vector{String}, table_rows::Vector{Vector{String}}, summary_items::Vector{Pair{String,String}})
    io = IOBuffer()
    println(io, "<!DOCTYPE html>")
    println(io, "<html lang=\"en\"><head><meta charset=\"utf-8\"><title>", escape_html(title), "</title>")
    println(io, "<style>body{font-family:-apple-system,BlinkMacSystemFont,\"Segoe UI\",sans-serif;margin:2rem;line-height:1.4}table{border-collapse:collapse;width:100%}th,td{border:1px solid #d0d7de;padding:0.5rem;text-align:left}th{background:#f6f8fa}code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}.meta{margin-bottom:1rem} .meta li{margin:0.25rem 0}</style>")
    println(io, "</head><body><h1>", escape_html(title), "</h1><ul class=\"meta\">")
    for (key, value) in summary_items
        println(io, "<li><strong>", escape_html(key), ":</strong> ", escape_html(value), "</li>")
    end
    println(io, "</ul><table><thead><tr>")
    for header in table_headers
        println(io, "<th>", escape_html(header), "</th>")
    end
    println(io, "</tr></thead><tbody>")
    for row in table_rows
        println(io, "<tr>")
        for cell in row
            println(io, "<td>", cell, "</td>")
        end
        println(io, "</tr>")
    end
    println(io, "</tbody></table></body></html>")
    String(take!(io))
end

function combined_summary_html(rows::Vector{Dict{String,Any}}, expected_shards::Vector{String}, completed_shards::Vector{String}, profile::AbstractString)
    missing = sort(setdiff(expected_shards, completed_shards))
    table_rows = [
        [
            "<code>$(escape_html(row["benchmark_id"]))</code>",
            "<code>$(escape_html(row["shard"]))</code>",
            escape_html(format_seconds_ns(row["time_ns_median"])),
            escape_html(format_bytes(row["memory_bytes"])),
            escape_html(string(row["allocs"])),
            "<code>$(escape_html(row["runner_kind"]))</code>",
            "<code>$(escape_html(row["hostname"]))</code>",
        ]
        for row in rows
    ]
    report_html(
        "Benchmark report",
        ["Benchmark ID", "Shard", "Median time", "Memory", "Allocs", "Runner", "Host"],
        table_rows,
        [
            "Profile" => profile,
            "Expected shards" => join(expected_shards, ", "),
            "Completed shards" => join(completed_shards, ", "),
            "Missing shards" => isempty(missing) ? "none" : join(missing, ", "),
        ],
    )
end

function combined_comparison_html(rows::Vector{Dict{String,Any}}, expected_shards::Vector{String}, completed_shards::Vector{String}, profile::AbstractString)
    missing = sort(setdiff(expected_shards, completed_shards))
    table_rows = [
        [
            "<code>$(escape_html(row["benchmark_id"]))</code>",
            "<code>$(escape_html(row["shard"]))</code>",
            escape_html(format_seconds_ns(row["target_time_ns_median"])),
            escape_html(format_seconds_ns(row["baseline_time_ns_median"])),
            escape_html(format_ratio(row["time_ratio"])),
            "<code>$(escape_html(row["judgement"]))</code>",
        ]
        for row in rows
    ]
    report_html(
        "Benchmark comparison report",
        ["Benchmark ID", "Shard", "Target time", "Baseline time", "Time ratio", "Judgement"],
        table_rows,
        [
            "Profile" => profile,
            "Expected shards" => join(expected_shards, ", "),
            "Completed shards" => join(completed_shards, ", "),
            "Missing shards" => isempty(missing) ? "none" : join(missing, ", "),
        ],
    )
end

function summary_files(input_dir::AbstractString, filename::AbstractString)
    matches = String[]
    for (root, _, files) in walkdir(input_dir)
        if filename in files
            push!(matches, joinpath(root, filename))
        end
    end
    sort!(matches)
end

function discovered_profile(rows::Vector{Dict{String,Any}}, default::AbstractString)
    isempty(rows) ? String(default) : String(rows[1]["profile"])
end

function shard_names(rows::Vector{Dict{String,Any}})
    sort!(unique!(String[row["shard"] for row in rows]))
end

function maybe_copy(path::AbstractString, env_name::AbstractString)
    destination = get(ENV, env_name, "")
    isempty(destination) && return nothing
    mkpath(dirname(destination))
    cp(path, destination; force=true)
    destination
end

function render_summary_report!(input_dir::AbstractString, output_dir::AbstractString; expected_shards::Vector{String}=String[])
    files = summary_files(input_dir, "summary.json")
    isempty(files) && error("No summary.json files were found under $input_dir")
    rows = Dict{String,Any}[]
    for file in files
        payload = JSON.parsefile(file)
        append!(rows, Dict{String,Any}.(payload["rows"]))
    end
    sort!(rows; by=row -> (row["benchmark_id"], row["shard"]))
    completed = shard_names(rows)
    expected = isempty(expected_shards) ? completed : expected_shards
    profile = discovered_profile(rows, get(ENV, "BENCHMARK_PROFILE", "small"))
    markdown = combined_summary_markdown(rows, expected, completed, profile)
    html = combined_summary_html(rows, expected, completed, profile)
    markdown_path = joinpath(output_dir, "benchmark-report.md")
    html_path = joinpath(output_dir, "benchmark-report.html")
    write(markdown_path, markdown)
    write(html_path, html)
    maybe_copy(markdown_path, "BENCHMARK_DOCS_MARKDOWN")
    maybe_copy(html_path, "BENCHMARK_DOCS_HTML")
    markdown_path, html_path
end

function render_comparison_report!(input_dir::AbstractString, output_dir::AbstractString; expected_shards::Vector{String}=String[])
    files = summary_files(input_dir, "comparison_summary.json")
    isempty(files) && error("No comparison_summary.json files were found under $input_dir")
    rows = Dict{String,Any}[]
    for file in files
        payload = JSON.parsefile(file)
        append!(rows, Dict{String,Any}.(payload["rows"]))
    end
    sort!(rows; by=row -> (row["benchmark_id"], row["shard"]))
    completed = shard_names(rows)
    expected = isempty(expected_shards) ? completed : expected_shards
    profile = discovered_profile(rows, get(ENV, "BENCHMARK_PROFILE", "small"))
    markdown = combined_comparison_markdown(rows, expected, completed, profile)
    html = combined_comparison_html(rows, expected, completed, profile)
    markdown_path = joinpath(output_dir, "benchmark-comparison-report.md")
    html_path = joinpath(output_dir, "benchmark-comparison-report.html")
    write(markdown_path, markdown)
    write(html_path, html)
    maybe_copy(markdown_path, "BENCHMARK_DOCS_MARKDOWN")
    maybe_copy(html_path, "BENCHMARK_DOCS_HTML")
    markdown_path, html_path
end

function run_benchmarks_from_env!()
    profile = normalize_profile(get(ENV, "BENCHMARK_PROFILE", "small"))
    shard = get(ENV, "BENCHMARK_SHARD", "all")
    result_dir = ensure_dir(get(ENV, "BENCHMARK_RESULT_DIR", default_result_dir()))
    shards = selected_shards(profile, shard)
    seconds = parse_float_env("BENCHMARK_SECONDS", 5.0)
    samples = parse_int_env("BENCHMARK_SAMPLES", 25)
    for shard_name in shards
        shard_dir = ensure_dir(joinpath(result_dir, shard_name))
        suite = filtered_suite(profile, shard_name)
        metadata = benchmark_metadata(profile, shard_name)
        results = BenchmarkTools.run(suite, verbose=true; seconds=seconds, samples=samples)
        BenchmarkTools.save(joinpath(shard_dir, "suite.json"), results)
        rows = summary_rows(results, metadata)
        write_json(joinpath(shard_dir, "metadata.json"), metadata)
        write_json(joinpath(shard_dir, "summary.json"), Dict("metadata" => metadata, "rows" => rows))
        write(joinpath(shard_dir, "summary.md"), summary_markdown(rows, metadata))
    end
    if parse_bool_env("BENCHMARK_RENDER_REPORT", true)
        render_summary_report!(result_dir, result_dir; expected_shards=shards)
    end
    println("Saved benchmark artifacts under $(abspath(result_dir))")
    nothing
end

function compare_benchmarks_from_env!()
    profile = normalize_profile(get(ENV, "BENCHMARK_PROFILE", "small"))
    shard = get(ENV, "BENCHMARK_SHARD", "all")
    result_dir = ensure_dir(get(ENV, "BENCHMARK_RESULT_DIR", default_result_dir()))
    baseline_ref = get(ENV, "BENCHMARK_BASELINE_REF", "")
    isempty(baseline_ref) && error("BENCHMARK_BASELINE_REF must be set for compare mode.")
    target_ref = get(ENV, "BENCHMARK_TARGET_REF", "")
    target_id = isempty(target_ref) ? nothing : target_ref
    shards = selected_shards(profile, shard)
    julia_cmd = `$(joinpath(Sys.BINDIR, Base.julia_exename())) --project=bench`
    for shard_name in shards
        shard_dir = ensure_dir(joinpath(result_dir, shard_name))
        bench_env = Dict(
            "BENCHMARK_PROFILE" => profile,
            "BENCHMARK_SHARD" => shard_name,
        )
        target_cfg = BenchmarkConfig(id=target_id, juliacmd=julia_cmd, env=bench_env)
        baseline_cfg = BenchmarkConfig(id=baseline_ref, juliacmd=julia_cmd, env=bench_env)
        target_path = joinpath(shard_dir, "target.json")
        baseline_path = joinpath(shard_dir, "baseline.json")
        target_results = benchmarkpkg(repo_root(), target_cfg; resultfile=target_path, verbose=true)
        baseline_results = benchmarkpkg(repo_root(), baseline_cfg; resultfile=baseline_path, verbose=true)
        judgement = judge(target_results, baseline_results, median; judgekwargs=Dict(:time_tolerance => 0.05))
        comparison_metadata = benchmark_metadata(profile, shard_name)
        comparison_metadata["target_ref"] = target_id === nothing ? "working-tree" : String(target_id)
        comparison_metadata["baseline_ref"] = baseline_ref
        rows = comparison_rows(
            PkgBenchmark.benchmarkgroup(target_results),
            PkgBenchmark.benchmarkgroup(baseline_results),
            PkgBenchmark.benchmarkgroup(judgement),
            comparison_metadata,
        )
        export_markdown(joinpath(shard_dir, "comparison.md"), judgement; export_invariants=true)
        write_json(joinpath(shard_dir, "metadata.json"), comparison_metadata)
        write_json(joinpath(shard_dir, "comparison_summary.json"), Dict("metadata" => comparison_metadata, "rows" => rows))
        write(joinpath(shard_dir, "comparison_summary.md"), comparison_markdown(rows, comparison_metadata))
    end
    if parse_bool_env("BENCHMARK_RENDER_REPORT", true)
        render_comparison_report!(result_dir, result_dir; expected_shards=shards)
    end
    println("Saved comparison artifacts under $(abspath(result_dir))")
    nothing
end

function render_report_from_env!()
    input_dir = get(ENV, "BENCHMARK_INPUT_DIR", get(ENV, "BENCHMARK_RESULT_DIR", default_result_dir()))
    output_dir = ensure_dir(get(ENV, "BENCHMARK_OUTPUT_DIR", input_dir))
    mode = get(ENV, "BENCHMARK_REPORT_MODE", "auto")
    expected = let raw = get(ENV, "BENCHMARK_EXPECTED_SHARDS", "")
        isempty(raw) ? String[] : split(raw, ',')
    end
    if mode == "comparison" || (mode == "auto" && !isempty(summary_files(input_dir, "comparison_summary.json")))
        render_comparison_report!(input_dir, output_dir; expected_shards=expected)
    else
        render_summary_report!(input_dir, output_dir; expected_shards=expected)
    end
    println("Rendered benchmark report under $(abspath(output_dir))")
    nothing
end

function write_placeholder_report!(path::AbstractString)
    mkpath(dirname(path))
    write(path, """
# Benchmark report

This page is generated from benchmark artifacts and intentionally ships with a placeholder in the repository.

To refresh it locally, run:

```julia
julia --project=bench bench/run_benchmarks.jl
BENCHMARK_INPUT_DIR=bench/results BENCHMARK_OUTPUT_DIR=bench/results BENCHMARK_DOCS_MARKDOWN=docs/src/benchmark_report.md julia --project=bench bench/render_report.jl
```
""")
    path
end
