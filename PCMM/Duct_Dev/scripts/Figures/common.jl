# Shared setup for every FigureN/simulate.jl (include at the top). Provides:
#   inputs                              the Duct_Dev_Figs project; each figure runs the DEFAULT model and
#                                       only lists the parameters it changes
#   save_snapshot(sim_id, sel, dst)     copy a snapshot picked by sel = :initial, :final, or a cell count
#   save_frame(sim_id, idx, dst)        copy snapshot number idx (or :final); frame_at(minutes) converts a time
#   fresh_dir(path)                     empty a figure's Outputs folder (keeps .gitkeep)
using PhysiCellModelManager, Printf

# Each figure runs as its own fresh `julia` process (no persisted project state), so initialize here.
initializeModelManager(joinpath(@__DIR__, "..", ".."))            # -> Duct_Dev
include(joinpath(@__DIR__, "..", "ParameterOptimization.jl"))     # -> output_index_at_cell_count

config_folder = custom_code_folder = rulesets_collection_folder = "Duct_Dev_Figs"
inputs = InputFolders(config_folder, custom_code_folder; rulesets_collection = rulesets_collection_folder)

# SVG (and full_data) save interval in the Duct_Dev_Figs config — frame index = minutes ÷ SVG_INTERVAL
const SVG_INTERVAL = 6
frame_at(minutes) = minutes ÷ SVG_INTERVAL

snapshot_path(sim_id, idx) = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output",
                                      idx === :final ? "final.svg" : @sprintf("snapshot%08d.svg", idx))

# Empty (or create) a figure's Outputs folder so stale snapshots never end up in a figure (keeps .gitkeep)
function fresh_dir(path)
    mkpath(path)
    for f in readdir(path)
        f == ".gitkeep" || rm(joinpath(path, f); force = true, recursive = true)
    end
    return path
end

# Snapshot index for `sel`: :initial, :final, or a cell count (first frame reaching it).
# NOTE: cell counts come from full_data, so this assumes SVG and full_data save intervals match (both 6 min).
function snapshot_index(sim_id, sel)
    sel === :initial && return 0
    sel === :final && return :final
    res = output_index_at_cell_count(sim_id, sel)
    if res === nothing
        @warn "sim $sim_id never reached $sel cells — using final frame"
        return :final
    end
    return first(res)
end

# Copy snapshot number `idx` (or :final) to `dst`
function save_frame(sim_id, idx, dst)
    cp(snapshot_path(sim_id, idx), dst; force = true)
    println("  -> ", dst)
end

# Copy the snapshot picked by `sel` to `dst`
save_snapshot(sim_id, sel, dst) = save_frame(sim_id, snapshot_index(sim_id, sel), dst)
