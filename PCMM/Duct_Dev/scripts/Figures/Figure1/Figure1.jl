using PhysiCellModelManager, Printf

# Fresh `julia` process has no persisted project state -- must initialize explicitly.
initializeModelManager(joinpath(@__DIR__, "..", ".."))

config_folder = custom_code_folder = rulesets_collection_folder = "Duct_Dev_Figs"
inputs = InputFolders(config_folder, custom_code_folder; rulesets_collection = rulesets_collection_folder)
resetDatabase(; force_reset=true, force_continue=true)



# Generating Snapshot for Figure 1A
dvs = [
    DiscreteVariation(configPath("max_time"), 7200),
    DiscreteVariation(configPath("CAF", "cycle", "rate", 0), 0.0),
    DiscreteVariation(configPath("Epithelial", "cycle", "rate", 0), 0.0),
]

out = run(inputs, dvs...; n_replicates = 1, force_recompile = true)
sim_id = only(simulationIDs(out))
println("sim_id = ", sim_id)

outdir = mkpath(joinpath(@__DIR__, "figure1"))
snapshot_path(idx) = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output",
                               idx === :final ? "final.svg" : @sprintf("snapshot%08d.svg", idx))

cp(snapshot_path(0),      joinpath(outdir, "A_initial.svg"); force = true)
cp(snapshot_path(:final), joinpath(outdir, "A_final.svg");   force = true)
println("  -> ", joinpath(outdir, "A_initial.svg"))
println("  -> ", joinpath(outdir, "A_final.svg"))

# Generating Snaptshot for Figure 1B

using PhysiCellModelManager, Printf

config_folder = custom_code_folder = rulesets_collection_folder = "Duct_Dev_Figs"
inputs = InputFolders(config_folder, custom_code_folder; rulesets_collection = rulesets_collection_folder)

dvs = [
    DiscreteVariation(configPath("max_time"), 7200),
    DiscreteVariation(configPath("user_parameters", "number_EP_cells"), 60),  # loose start — set to whatever you used before
    DiscreteVariation(configPath("Epithelial", "cycle", "rate", 0), 0.00072),
    DiscreteVariation(configPath("CAF", "cycle", "rate", 0), 0),  
]

out = run(inputs, dvs...; n_replicates = 1, force_recompile = true)
sim_id = only(simulationIDs(out))
println("sim_id = ", sim_id)

outdir = mkpath(joinpath(@__DIR__, "figure1"))
snapshot_path(idx) = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output",
                               idx === :final ? "final.svg" : @sprintf("snapshot%08d.svg", idx))

cp(snapshot_path(0),      joinpath(outdir, "B_initial.svg"); force = true)
cp(snapshot_path(:final), joinpath(outdir, "B_final.svg");   force = true)
println("  -> ", joinpath(outdir, "B_initial.svg"))
println("  -> ", joinpath(outdir, "B_final.svg"))

