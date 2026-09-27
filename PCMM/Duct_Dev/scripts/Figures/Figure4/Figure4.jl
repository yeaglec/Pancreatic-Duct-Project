using PhysiCellModelManager, Printf

# Fresh `julia` process has no persisted project state (unlike an interactive
# session where you'd already have called this once) -- must initialize explicitly.
initializeModelManager(joinpath(@__DIR__, "..", ".."))
include(joinpath(@__DIR__, "..", "ParameterOptimization.jl"))   # -> output_index_at_cell_count

config_folder = custom_code_folder = rulesets_collection_folder = "Duct_Dev_Figs"
inputs = InputFolders(config_folder, custom_code_folder; rulesets_collection = rulesets_collection_folder)

# Same fixed baseline as Figure 3's seeded-CAF proliferation run
BASELINE = [
    DiscreteVariation(configPath("max_time"), 10080),
    DiscreteVariation(configPath("user_parameters", "is_restore_lin"), 1),
    DiscreteVariation(configPath("user_parameters", "number_CAF_cells"), 1),
    DiscreteVariation(configPath("Epithelial", "cycle", "rate", 0), 0.0),
]

TARGET_CELL_COUNT = 175   # compare all variants at the same growth milestone

# (label, xml path, values to sweep) — baseline value included in each list as a reference point
SWEEPS = [
    ("kernel_sigma", configPath("user_parameters", "membrane_force_smoothing_sigma"), [0.0, 10.0, 25.0, 50.0]),
    ("strain",       configPath("user_parameters", "seg_lin"),                        [0.05, 0.1, 0.2, 0.4]),
    ("restoring",    configPath("user_parameters", "home_lin"),                       [0.0, 0.005, 0.01, 0.02]),
]

outdir_base = mkpath(joinpath(@__DIR__, "figure4"))
snapshot_path(sim_id, idx) = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output",
                                       idx === :final ? "final.svg" : @sprintf("snapshot%08d.svg", idx))

for (name, path, values) in SWEEPS
    outdir = mkpath(joinpath(outdir_base, name))
    for v in values
        dvs = [BASELINE..., DiscreteVariation(path, v)]
        out = run(inputs, dvs...; n_replicates = 1, force_recompile = false)
        sim_id = only(simulationIDs(out))

        res = output_index_at_cell_count(sim_id, TARGET_CELL_COUNT)
        idx = res === nothing ? :final : first(res)
        res === nothing && @warn "$name = $v never reached $TARGET_CELL_COUNT cells — using final frame"

        dst = joinpath(outdir, "v$(v).svg")
        cp(snapshot_path(sim_id, idx), dst; force = true)
        println(name, " = ", v, " -> sim_id ", sim_id, "  -> ", dst)
    end
end