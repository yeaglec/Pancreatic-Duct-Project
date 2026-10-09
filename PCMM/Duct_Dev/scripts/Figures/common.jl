# Shared setup for every FigureN/simulate.jl (include at the top). Provides:
#   inputs                              the Duct_Dev_Figs project; each figure runs the DEFAULT model and
#                                       only lists the parameters it changes
#   save_snapshot(sim_id, sel, dst)     copy a snapshot picked by sel = :initial, :final, or a cell count
#   save_frame(sim_id, idx, dst)        copy snapshot number idx (or :final); frame_at(minutes) converts a time
#   fresh_dir(path)                     empty a figure's Outputs folder (keeps .gitkeep)
#   frame_metrics(sim_id, frame)        the four duct metrics (Metrics.jl) of one saved frame; load_membrane / load_cells
using PhysiCellModelManager, Printf, DelimitedFiles

# Each figure runs as its own fresh `julia` process (no persisted project state), so initialize here.
initializeModelManager(joinpath(@__DIR__, "..", ".."))            # -> Duct_Dev
include(joinpath(@__DIR__, "..", "ParameterOptimization.jl"))     # -> output_index_at_cell_count
include(joinpath(@__DIR__, "..", "Metrics.jl"))

config_folder = custom_code_folder = rulesets_collection_folder = "Duct_Dev_Figs"
inputs = InputFolders(config_folder, custom_code_folder; rulesets_collection = rulesets_collection_folder)

# SVG (and full_data) save interval in the Duct_Dev_Figs config — frame index = minutes ÷ SVG_INTERVAL
const SVG_INTERVAL = 6
frame_at(minutes) = minutes ÷ SVG_INTERVAL

output_dir(sim_id) = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output")
snapshot_path(sim_id, idx) = joinpath(output_dir(sim_id), idx === :final ? "final.svg" : @sprintf("snapshot%08d.svg", idx))

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

# Membrane nodes (N×2) and live duct cells (x, y, radius, type) of saved frame `frame` (same 6-min index as the snapshots)
load_membrane(sim_id, frame) = readdlm(joinpath(output_dir(sim_id), "boundary_t$(frame).csv"), ',', Float64)
function load_cells(sim_id, frame)
    c = PhysiCellSnapshot(sim_id, frame; include_cells = true).cells
    c = c[(.!c.dead) .& (c.cell_type_name .!= "BM"), :]   # the basement-membrane agents are not part of the duct
    return c.position_1, c.position_2, c.radius, c.cell_type_name
end

# The four metrics (Metrics.jl) of one saved frame, plus the cell count
function frame_metrics(sim_id, frame; reach = 10.0)
    m = load_membrane(sim_id, frame)
    x, y, r, _ = load_cells(sim_id, frame)
    return (IC = inverse_circularity(m), indent = max_indent_depth(m),
            off_membrane = off_membrane_fraction(m, x, y, r; reach), filling = lumen_filling(m, r), cells = length(x))
end
