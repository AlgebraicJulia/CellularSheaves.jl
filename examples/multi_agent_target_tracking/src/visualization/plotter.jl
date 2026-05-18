module Plotter

using CSV
using DataFrames
using Plots
using Statistics
using JSON3
using Printf

const STATE_SUFFIX = "_state_data.csv"

# -----------------------------------------------------------------------------
# File/data helpers
# -----------------------------------------------------------------------------

function _extract_num(s::AbstractString)
    m = match(r"(\d+)\D*$", s)
    return m === nothing ? typemax(Int) : parse(Int, m.captures[1])
end

function _list_state_files(dir::AbstractString)
    isdir(dir) || return String[]
    return sort(filter(f -> endswith(f, STATE_SUFFIX), readdir(dir)); by=_extract_num)
end

_name_from_file(f::AbstractString) = replace(f, STATE_SUFFIX => "")
_read_state(path::AbstractString) = CSV.read(path, DataFrame)
_rms(x::AbstractVector{<:Real}) = sqrt(mean(abs2, x))

function _example_root()
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function _default_outdir()
    return joinpath(_example_root(), "simulation_data")
end

function _default_config_path()
    return joinpath(_example_root(), "configurations", "config_common.json")
end

function get_simulation_data(outdir::AbstractString)
    agent_dir = joinpath(outdir, "agent_data")
    target_dir = joinpath(outdir, "target_data")

    agent_files = _list_state_files(agent_dir)
    target_files = _list_state_files(target_dir)

    agent_names = _name_from_file.(agent_files)
    target_names = _name_from_file.(target_files)
    agent_dfs = [_read_state(joinpath(agent_dir, f)) for f in agent_files]
    target_dfs = [_read_state(joinpath(target_dir, f)) for f in target_files]

    return agent_names, agent_dfs, target_names, target_dfs
end

function _load_config(config_path::Union{Nothing,AbstractString})
    config_path === nothing && return nothing
    isfile(config_path) || return nothing
    return Dict{String,Any}(JSON3.read(read(config_path, String), Dict{String,Any}))
end

function _nested_to_matrix(x)
    rows = collect(x)
    isempty(rows) && return zeros(Float64, 0, 0)
    nrows = length(rows)
    ncols = length(rows[1])
    M = zeros(Float64, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        M[i, j] = Float64(rows[i][j])
    end
    return M
end

function _edge_pairs(x)
    pairs = Tuple{Int,Int}[]
    for e in collect(x)
        length(e) >= 2 || continue
        push!(pairs, (Int(e[1]), Int(e[2])))
    end
    return pairs
end

function _agent_neighbors(n_agents::Int, edge_set)
    neighbors = [Int[] for _ in 1:n_agents]
    for (i, j) in _edge_pairs(edge_set)
        if 1 <= i <= n_agents && 1 <= j <= n_agents
            push!(neighbors[i], j)
            push!(neighbors[j], i)
        end
    end
    return neighbors
end

function _chunk_triplets(v::Vector{Int})
    groups = Vector{Int}[]
    i = 1
    while i + 2 <= length(v)
        push!(groups, v[i:i+2])
        i += 3
    end
    return groups
end

function _pos_at(df::DataFrame, row::Int)
    r = clamp(row, 1, nrow(df))
    return (Float64(df[r, Symbol("Position X")]), Float64(df[r, Symbol("Position Y")]))
end

function _all_xy_at(agent_dfs, target_dfs, rows::Vector{Int})
    pts = Tuple{Float64,Float64}[]
    for r in rows
        for df in agent_dfs
            push!(pts, _pos_at(df, r))
        end
        for df in target_dfs
            push!(pts, _pos_at(df, r))
        end
    end
    return pts
end

function _xy_limits(agent_dfs, target_dfs, rows::Vector{Int})
    pts = _all_xy_at(agent_dfs, target_dfs, rows)
    isempty(pts) && return (-1.0, 1.0), (-1.0, 1.0)

    xs = first.(pts)
    ys = last.(pts)
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)

    # Use a square span so both axes have the same range, then pad.
    # This replaces aspect_ratio=:equal (which fights explicit xlim/ylim
    # inside subplot layouts in Plots.jl) with manually enforced square limits.
    cx = (xmin + xmax) / 2
    cy = (ymin + ymax) / 2
    half = max(xmax - xmin, ymax - ymin, 1e-9) / 2
    pad = 0.12 * half * 2
    lo = -half - pad
    hi =  half + pad
    return (cx + lo, cx + hi), (cy + lo, cy + hi)
end

function _draw_line!(p, subplot_id::Int, a, b; label="", kwargs...)
    plot!(p, [a[1], b[1]], [a[2], b[2]]; subplot=subplot_id, label=label, kwargs...)
end

# -----------------------------------------------------------------------------
# Plot 1: synchronization/tracking error norm
# -----------------------------------------------------------------------------

function plot_error_norm(agent_names, agent_dfs; show::Bool=true, save_path=nothing)
    isempty(agent_dfs) && return nothing

    p = plot(; xlabel="Time (s)",
               ylabel="Synchronization Error Norm (m)",
               title="Tracking error norm",
               yscale=:log10,
               legend=:topright,
               legendfontsize=7,
               size=(800, 480))

    rows = collect(zip(agent_names, agent_dfs))
    sort!(rows; by=x -> _extract_num(x[1]))

    all_vals = Float64[]
    for (name, df) in rows
        err = Float64.(df[!, Symbol("Synchronization Error Norm")])
        append!(all_vals, err)
        rms_value = _rms(err)
        agent_num = _extract_num(name)
        label = agent_num == typemax(Int) ?
            "$(name): RMS $(round(rms_value, digits=2)) m" :
            "Agent $(agent_num): RMS $(round(rms_value, digits=2)) m"
        plot!(p, df[!, :Time], err; label=label, lw=1.5)
    end

    # Clamp y-axis so a near-zero transient doesn't collapse the scale.
    pos_vals = filter(>(0), all_vals)
    if !isempty(pos_vals)
        yfloor = max(1e-4, 10^floor(log10(minimum(pos_vals))))
        yceil  = 10^ceil(log10(maximum(all_vals))) * 2.0
        ylims!(p, yfloor, yceil)
    end

    if !isempty(agent_dfs)
        tmax = maximum(agent_dfs[1][!, :Time])
        xlims!(p, 0.0, min(10.0, Float64(tmax)))
    end

    save_path !== nothing && savefig(p, save_path)
    show && display(p)
    return p
end

# -----------------------------------------------------------------------------
# Plot 2: 3D trajectories
# -----------------------------------------------------------------------------

function plot_trajectories_3d(agent_names, agent_dfs, target_names, target_dfs;
                              show::Bool=true, save_path=nothing)

    p = plot(
        xlabel="X (m)",
        ylabel="Y (m)",
        zlabel="Z (m)",
        title="Trajectories",

        # Bigger figure
        size=(1100, 800),

        # Smaller legend
        legend=:outerleft,
        legendfontsize=6,

        # Better spacing
        left_margin=8Plots.mm,
        right_margin=12Plots.mm,
        bottom_margin=12Plots.mm,
        top_margin=8Plots.mm,

        # Better 3D viewing angle
        camera=(45, 25)
    )

    for (name, df) in zip(agent_names, agent_dfs)
        plot!(p,
              df[!, Symbol("Position X")],
              df[!, Symbol("Position Y")],
              df[!, Symbol("Position Z")];
              label=name, lw=2, linestyle=:solid)
    end

    for (name, df) in zip(target_names, target_dfs)
        plot!(p,
              df[!, Symbol("Position X")],
              df[!, Symbol("Position Y")],
              df[!, Symbol("Position Z")];
              label=name, lw=3, linestyle=:dash)
    end

    save_path !== nothing && savefig(p, save_path)
    show && display(p)
    return p
end

# -----------------------------------------------------------------------------
# Plot 3: 2D snapshots, initial and final, with TT / AT / AA links
# -----------------------------------------------------------------------------

function plot_snapshots_2d(agent_names, agent_dfs, target_names, target_dfs;
                           config_path::Union{Nothing,AbstractString}=nothing,
                           show::Bool=true,
                           save_path=nothing)
    isempty(agent_dfs) && isempty(target_dfs) && return nothing

    cfg = _load_config(config_path)
    n_agents = length(agent_dfs)
    n_targets = length(target_dfs)

    target_edge_set = cfg === nothing ? Tuple{Int,Int}[] : _edge_pairs(get(cfg, "target_edge_set", Any[]))
    agent_edge_set_raw = cfg === nothing ? Any[] : get(cfg, "agent_edge_set", Any[])
    agent_neighbors = _agent_neighbors(n_agents, agent_edge_set_raw)
    pinning_matrix = cfg === nothing || !haskey(cfg, "pinning_matrix") ? zeros(Float64, 0, 0) : _nested_to_matrix(cfg["pinning_matrix"])

    k_square = 4
    square_group = n_agents >= k_square ? collect(n_agents-k_square+1:n_agents) : Int[]
    tri_pool = n_agents >= k_square ? collect(1:n_agents-k_square) : collect(1:n_agents)
    tri_groups = _chunk_triplets(tri_pool)

    max_len = maximum([nrow(df) for df in vcat(agent_dfs, target_dfs)])
    rows = [1, max_len]
    titles = ["t = 0", "t = final"]

    # Same axis limits for both plots, so initial/final are comparable.
    xlim, ylim = _xy_limits(agent_dfs, target_dfs, rows)

    tri_colors = [:purple, :green, :orange, :brown]

    function draw_single_snapshot(row::Int, title_text::String; legend_on::Bool=true)
        p = plot(
        xlabel="X (m)",
        ylabel="Y (m)",
        title=title_text,
        grid=true,
        legend=(legend_on ? :outerright : false),
        legendfontsize=6,
        legendtitle="Formations",
        legendtitlefontsize=7,
        size=(850, 650),
        left_margin=10Plots.mm,
        right_margin=35Plots.mm,
        bottom_margin=14Plots.mm,
        top_margin=8Plots.mm
        )

        did_tt = false
        did_at = false
        did_square = false
        did_inter = false
        did_tri = falses(length(tri_groups))

        # Target-target links
        for (i, j) in target_edge_set
            if 1 <= i <= n_targets && 1 <= j <= n_targets
                a = _pos_at(target_dfs[i], row)
                b = _pos_at(target_dfs[j], row)

                plot!(p, [a[1], b[1]], [a[2], b[2]];
                    color=:red,
                    lw=2.2,
                    alpha=0.85,
                    label="")

                did_tt = true
            end
        end

        # Dummy legend entry for target formation
        if legend_on && did_tt
            plot!(p, [NaN], [NaN];
                color=:red,
                lw=2.2,
                alpha=0.85,
                label="Target Formation")
        end

        # Agent-target pinning links
        if size(pinning_matrix, 1) == n_agents
            for ai in 1:size(pinning_matrix, 1), tj in 1:size(pinning_matrix, 2)
                if pinning_matrix[ai, tj] != 0.0 && 1 <= tj <= n_targets
                    a = _pos_at(agent_dfs[ai], row)
                    b = _pos_at(target_dfs[tj], row)

                    plot!(p, [a[1], b[1]], [a[2], b[2]];
                          color=:cyan,
                          lw=1.5,
                          linestyle=:dot,
                          alpha=0.80,
                        #   label=(legend_on && !did_at) ? "Agent-Target" : "")
                          label="")

                    did_at = true
                end
            end
        end

        # Agent-agent links inside triangle groups
        for (gidx, g) in enumerate(tri_groups)
            col = tri_colors[mod1(gidx, length(tri_colors))]

            for i in g
                for j in agent_neighbors[i]
                    if j > i && j in g
                        a = _pos_at(agent_dfs[i], row)
                        b = _pos_at(agent_dfs[j], row)

                        plot!(p, [a[1], b[1]], [a[2], b[2]];
                              color=col,
                              lw=1.8,
                              alpha=0.85,
                              label=(legend_on && !did_tri[gidx]) ? "Triangle Formation $(gidx)" : "")

                        did_tri[gidx] = true
                    end
                end
            end
        end

        # Agent-agent links inside square group
        for i in square_group
            for j in agent_neighbors[i]
                if j > i && j in square_group
                    a = _pos_at(agent_dfs[i], row)
                    b = _pos_at(agent_dfs[j], row)

                    plot!(p, [a[1], b[1]], [a[2], b[2]];
                          color=:blue,
                          lw=1.8,
                          alpha=0.85,
                          label=(legend_on && !did_square) ? "Square Formation" : "")

                    did_square = true
                end
            end
        end

        # Agent-agent links across formations
        labels_map = Dict{Int,String}()

        for i in square_group
            labels_map[i] = "square"
        end

        for (gidx, g) in enumerate(tri_groups)
            for i in g
                labels_map[i] = "tri$(gidx)"
            end
        end

        for i in 1:n_agents
            for j in agent_neighbors[i]
                if j > i && haskey(labels_map, i) && haskey(labels_map, j) && labels_map[i] != labels_map[j]
                    a = _pos_at(agent_dfs[i], row)
                    b = _pos_at(agent_dfs[j], row)

                    plot!(p, [a[1], b[1]], [a[2], b[2]];
                          color=:olive,
                          lw=1.5,
                          linestyle=:dot,
                          alpha=0.80,
                        #   label=(legend_on && !did_inter) ? "Agent-Agent (out of formation)" : "")
                            label="")

                    did_inter = true
                end
            end
        end

        # Targets
        for (name, df) in zip(target_names, target_dfs)
            x, y = _pos_at(df, row)
            scatter!(p, [x], [y];
                     markershape=:x,
                     markersize=8,
                     markerstrokewidth=2,
                     label="")
        end

        # Agents
        for (name, df) in zip(agent_names, agent_dfs)
            x, y = _pos_at(df, row)
            scatter!(p, [x], [y];
                     markershape=:circle,
                     markersize=5,
                     markerstrokecolor=:black,
                     markerstrokewidth=0.5,
                     label="")
        end

        xlims!(p, xlim...)
        ylims!(p, ylim...)

        return p
    end

    p_initial = draw_single_snapshot(rows[1], titles[1]; legend_on=true)
    p_final = draw_single_snapshot(rows[2], titles[2]; legend_on=true)

    # p_initial = draw_single_snapshot(rows[1], titles[1]; legend_on=false)
    # p_final = draw_single_snapshot(rows[2], titles[2]; legend_on=false)

    if save_path !== nothing
        save_dir = dirname(save_path)
        mkpath(save_dir)

        initial_path = joinpath(save_dir, "snapshot_t0.png")
        final_path   = joinpath(save_dir, "snapshot_final.png")

        savefig(p_initial, initial_path)
        savefig(p_final, final_path)

        @info "Saved initial 2D snapshot" initial_path
        @info "Saved final 2D snapshot" final_path
    end

    if show
        display(p_initial)
        display(p_final)
    end

    return (p_initial=p_initial, p_final=p_final)
end

# -----------------------------------------------------------------------------
# Average error utilities
# -----------------------------------------------------------------------------

function compute_avg_error(outdir::AbstractString; save::Bool=true)
    agent_names, agent_dfs, _, _ = get_simulation_data(outdir)
    isempty(agent_dfs) && error("No agent CSV files found in $(joinpath(outdir, "agent_data"))")

    time = agent_dfs[1][!, :Time]
    avg = zeros(Float64, length(time))

    for df in agent_dfs
        err = df[!, Symbol("Synchronization Error Norm")]
        @assert length(err) == length(time) "All agent CSVs must have the same length."
        avg .+= Float64.(err)
    end
    avg ./= length(agent_dfs)

    result = DataFrame(Time=time, avg_error=avg)
    if save
        path = joinpath(outdir, "avg_error.csv")
        CSV.write(path, result)
        @info "Saved average error CSV" path
    end
    return result
end

function _label_from_run_dir(run_dir::AbstractString)
    base = basename(normpath(run_dir))
    m = match(r"k[_=]([0-9\.eE+-]+)", base)
    return m === nothing ? base : m.captures[1]
end

function plot_avg_errors(run_dirs::Vector{<:AbstractString};
                         show::Bool=true,
                         save_path::Union{Nothing,AbstractString}=nothing,
                         save_csv::Union{Nothing,AbstractString}=nothing)
    isempty(run_dirs) && error("Pass at least one run directory.")

    p = plot(; xlabel="Time (s)", ylabel="Average Tracking Error", yscale=:log10,
               legend=:best, legendfontsize=7)
    wide = nothing

    for run_dir in sort(String.(run_dirs))
        avg_path = joinpath(run_dir, "avg_error.csv")
        df = isfile(avg_path) ? CSV.read(avg_path, DataFrame) : compute_avg_error(run_dir; save=true)
        label = _label_from_run_dir(run_dir)
        plot!(p, df[!, :Time], df[!, :avg_error]; label="k=$(label)", lw=2)

        renamed = DataFrame(Time=df[!, :Time])
        renamed[!, Symbol("k=$(label)")] = df[!, :avg_error]
        wide = wide === nothing ? renamed : outerjoin(wide, renamed, on=:Time)
    end

    save_csv !== nothing && wide !== nothing && CSV.write(save_csv, sort(wide, :Time))
    save_path !== nothing && savefig(p, save_path)
    show && display(p)
    return p
end

function plot_avg_errors_from_parent(parent_dir::AbstractString; kwargs...)
    isdir(parent_dir) || error("Parent directory does not exist: $(parent_dir)")
    run_dirs = [joinpath(parent_dir, d) for d in readdir(parent_dir) if isdir(joinpath(parent_dir, d))]
    return plot_avg_errors(run_dirs; kwargs...)
end

# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

function plot_from_csv(outdir::AbstractString;
                       config_path::Union{Nothing,AbstractString}=nothing,
                       show::Bool=true,
                       save_dir=nothing,
                       snapshots::Bool=true)
    agent_names, agent_dfs, target_names, target_dfs = get_simulation_data(outdir)

    if isempty(agent_dfs) && isempty(target_dfs)
        @warn "No simulation data found." outdir
        return nothing
    end

    save_dir !== nothing && mkpath(save_dir)

    p_err = plot_error_norm(agent_names, agent_dfs;
                            show=show,
                            save_path=save_dir === nothing ? nothing : joinpath(save_dir, "error_norm.png"))

    p_traj = plot_trajectories_3d(agent_names, agent_dfs, target_names, target_dfs;
                                  show=show,
                                  save_path=save_dir === nothing ? nothing : joinpath(save_dir, "trajectories_3d.png"))

    p_snap = nothing
    if snapshots
        snap_path = save_dir === nothing ? nothing : joinpath(save_dir, "snapshots_2d.png")
        @info "Generating 2D snapshots" snap_path config_path n_agents=length(agent_dfs) n_targets=length(target_dfs)
        p_snap = plot_snapshots_2d(agent_names, agent_dfs, target_names, target_dfs;
                                   config_path=config_path,
                                   show=show,
                                   save_path=snap_path)
    else
        @info "Skipping 2D snapshots because snapshots=false"
    end

    return (p_err=p_err, p_traj=p_traj, p_snap=p_snap)
end

function results(; outdir=nothing,
                   config_path=nothing,
                   show::Bool=false,
                   save_dir=nothing,
                   snapshots::Bool=true)
    outdir === nothing && (outdir = _default_outdir())
    config_path === nothing && (config_path = _default_config_path())
    save_dir === nothing && (save_dir = joinpath(_example_root(), "plots"))

    return plot_from_csv(
        outdir;
        config_path=config_path,
        show=show,
        save_dir=save_dir,
        snapshots=snapshots
    )
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    using .Plotter
    Plotter.results()
end