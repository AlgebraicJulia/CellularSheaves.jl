using CellularSheaves
using Plots
using LaTeXStrings

const data_root   = joinpath(@__DIR__, "../data/asynch")
const figures_dir = joinpath(@__DIR__, "../figures/asynch")
mkpath(figures_dir)

subsample(v, n=2000) = v[1:max(1, length(v) ÷ n):end]

# ── convergence_vs_delay ──────────────────────────────────────────────────────

let dr = joinpath(data_root, "converge_vs_delay")
    Bs     = [1, 10, 50, 100, 200]
    labels = ["B=1", "B=10", "B=50", "B=100", "B=200"]
    topologies      = ["regular_graph", "er_graph", "complete_graph"]
    topology_titles = ["4-Regular", "Erdős–Rényi", "Star"]

    loss_plots = map(zip(topologies, topology_titles)) do (top, title)
        plt = empty_experiment_plot(; title=title)
        for (i, (_, lab)) in enumerate(zip(Bs, labels))
            path = joinpath(dr, top, "loss$i.csv")
            isfile(path) || continue
            plot_log_loss_curve!(plt, subsample(vec(load_trajectory(path))), lab)
        end
        plt
    end

    error_plots = map(zip(topologies, topology_titles)) do (top, _)
        plt = empty_experiment_plot(; xlabel=L"\mathrm{Iteration},\quad t")
        for (i, (_, lab)) in enumerate(zip(Bs, labels))
            path = joinpath(dr, top, "error$i.csv")
            isfile(path) || continue
            plot_log_loss_curve!(plt, subsample(vec(load_trajectory(path))), lab)
        end
        plt
    end

    l = @layout [a b c; d e f]
    plt = plot(
        loss_plots..., error_plots...,
        layout=l,
        size=(1900, 700),
        ylabel=[L"f(\mathbf{x})" "" "" L"\|\mathbf{x}-\mathbf{x}^*\|/\|\mathbf{x}^*\|" "" ""],
        legend=:topright,
    )
    savefig(plt, joinpath(figures_dir, "convergence_vs_delay.svg"))
    @info "Saved convergence_vs_delay.svg"
end

# ── step_size_comparison ──────────────────────────────────────────────────────

let dr = joinpath(data_root, "ss_vs_converge")
    plt = empty_experiment_plot(;
        xlabel=L"\mathrm{Iteration},\quad t",
        ylabel=L"f(\mathbf{x})",
        legend=:topright,
    )
    for (fname, lab) in [
        ("loss_synch.csv",     "Synch  1/K"),
        ("loss_synch_ss.csv",  "Asynch 1/K"),
        ("loss_asynch_ss.csv", "Asynch 0.1/K"),
    ]
        path = joinpath(dr, fname)
        isfile(path) || continue
        plot_log_loss_curve!(plt, subsample(vec(load_trajectory(path))), lab)
    end
    savefig(plt, joinpath(figures_dir, "step_size_comparison.svg"))
    @info "Saved step_size_comparison.svg"
end

# ── restriction_map_comparison ────────────────────────────────────────────────

let dr = joinpath(data_root, "rm_experiment")
    Bs     = [1, 10, 50, 100, 200]
    labels = ["B=1", "B=10", "B=50", "B=100", "B=200"]
    types       = ["constant_sheaf", "vector_bundle", "matrix_weighted"]
    type_titles = ["Constant Sheaf", "Random Vector Bundle", "Matrix-Weighted"]

    function load_curves(subdir, prefix)
        map(1:length(Bs)) do i
            path = joinpath(dr, subdir, "$(prefix)$i.csv")
            isfile(path) ? subsample(vec(load_trajectory(path))) : Float64[]
        end
    end

    loss_plots = map(zip(types, type_titles)) do (t, title)
        plt = empty_experiment_plot(; title=title)
        for (curves, lab) in zip(load_curves(t, "loss"), labels)
            isempty(curves) && continue
            plot_log_loss_curve!(plt, curves, lab)
        end
        plt
    end

    error_plots = map(types) do t
        plt = empty_experiment_plot()
        for (curves, lab) in zip(load_curves(t, "error"), labels)
            isempty(curves) && continue
            plot_log_loss_curve!(plt, curves, lab)
        end
        plt
    end

    proj_plots = map(types) do t
        plt = empty_experiment_plot(; xlabel=L"\mathrm{Iteration},\quad t")
        for (curves, lab) in zip(load_curves(t, "proj_error"), labels)
            isempty(curves) && continue
            plot_log_loss_curve!(plt, curves, lab)
        end
        plt
    end

    l = @layout [a b c; d e f; g h i]
    plt = plot(
        loss_plots..., error_plots..., proj_plots...,
        layout=l,
        size=(1900, 1000),
        ylabel=[L"f(\mathbf{x})" "" "" L"\|\mathbf{x}-\mathbf{x}^*\|/\|\mathbf{x}^*\|" "" "" L"\|\mathbf{x}-\Pi[\mathbf{x}(0)]\|/\|\Pi[\mathbf{x}(0)]\|" "" ""],
        legend=:topright,
    )
    savefig(plt, joinpath(figures_dir, "restriction_map_comparison.svg"))
    @info "Saved restriction_map_comparison.svg"
end

# ── orthogonal_projection ─────────────────────────────────────────────────────

let dr = joinpath(data_root, "op_experiment"),
    results_path = joinpath(data_root, "op_experiment", "results.csv")
    if isfile(results_path)
        results = vec(load_trajectory(results_path))
        Bs = [2^(i - 1) for i in 1:length(results)]
        plt = scatter(
            Bs, results,
            xscale=:log2,
            xlabel=L"B",
            ylabel=L"\|\mathbf{x}^* - \mathbf{x}(T)^\perp\|",
            label="",
            thickness_scaling=1.9,
            title="Convergence error vs. max delay",
        )
        savefig(plt, joinpath(figures_dir, "orthogonal_projection.svg"))
        @info "Saved orthogonal_projection.svg"
    else
        @warn "Data not found at $results_path. Run docs/scripts/acc26_experiments.jl first."
    end
end
