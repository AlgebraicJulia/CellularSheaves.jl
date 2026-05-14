# src/io/data_manager.jl

mutable struct DataManager
    outdir::String
    agent_ios::Dict{String, IO}
    target_ios::Dict{String, IO}
end

function DataManager(outdir::AbstractString="simulation_data")
    return DataManager(String(outdir), Dict{String,IO}(), Dict{String,IO}())
end

function _ensure_dirs(dm::DataManager)
    mkpath(joinpath(dm.outdir, "agent_data"))
    mkpath(joinpath(dm.outdir, "target_data"))
end

function init_data_manager!(dm::DataManager, agents, targets)
    _ensure_dirs(dm)

    for ag in agents
        id = String(ag.id)
        path = joinpath(dm.outdir, "agent_data", "$(id)_state_data.csv")
        io = open(path, "w")
        write(io, "Time,Position X,Position Y,Position Z,Synchronization Error Norm\n")
        dm.agent_ios[id] = io
    end

    for tg in targets
        id = String(tg.id)
        path = joinpath(dm.outdir, "target_data", "$(id)_state_data.csv")
        io = open(path, "w")
        write(io, "Time,Position X,Position Y,Position Z\n")
        dm.target_ios[id] = io
    end

    return dm
end

function _csv_float(x)
    return string(Float64(x))
end

function _write_agent_row!(io::IO, agent, col::Int, time::Float64)
    err_norm = sqrt(sum(abs2, agent.synchronization_error))
    write(io, _csv_float(time), ",",
              _csv_float(agent.positions[1, col]), ",",
              _csv_float(agent.positions[2, col]), ",",
              _csv_float(agent.positions[3, col]), ",",
              _csv_float(err_norm), "\n")
end

function _write_target_row!(io::IO, target, col::Int, time::Float64)
    write(io, _csv_float(time), ",",
              _csv_float(target.positions[1, col]), ",",
              _csv_float(target.positions[2, col]), ",",
              _csv_float(target.positions[3, col]), "\n")
end

"""
    save_step!(dm, agents, targets, step, dt; python_style=true)

Writes Python-compatible CSV files. With `python_style=true`, the saved row uses
state column `step - 1`, matching Python's `save_state_to_csv(step, ...)` behavior.
"""
function save_step!(dm::DataManager, agents, targets, step::Int, dt::Float64; python_style::Bool=true)
    col = python_style ? max(step - 1, 1) : step
    time = (step - 1) * dt

    for ag in agents
        io = dm.agent_ios[String(ag.id)]
        _write_agent_row!(io, ag, col, time)
    end
    for tg in targets
        io = dm.target_ios[String(tg.id)]
        _write_target_row!(io, tg, col, time)
    end
    return nothing
end

function close_all_files!(dm::DataManager)
    for (_, io) in dm.agent_ios
        close(io)
    end
    for (_, io) in dm.target_ios
        close(io)
    end
    empty!(dm.agent_ios)
    empty!(dm.target_ios)
    return nothing
end
