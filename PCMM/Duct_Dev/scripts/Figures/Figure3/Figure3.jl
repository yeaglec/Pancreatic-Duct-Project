using PhysiCellModelManager, Printf

# Fresh `julia` process has no persisted project state (unlike an interactive
# session where you'd already have called this once) -- must initialize explicitly.
initializeModelManager(joinpath(@__DIR__, "..", ".."))
include(joinpath(@__DIR__, "..", "ParameterOptimization.jl"))   # -> output_index_at_cell_count only; no report/HTML machinery

config_folder = custom_code_folder = rulesets_collection_folder = "Duct_Dev_Figs"
inputs = InputFolders(config_folder, custom_code_folder; rulesets_collection = rulesets_collection_folder)

dvs = [
    DiscreteVariation(configPath("max_time"), 10080),
    DiscreteVariation(configPath("user_parameters", "is_restore_lin"), 1),   # not redundant — base default is 0
    DiscreteVariation(configPath("user_parameters", "number_CAF_cells"), 1),
    DiscreteVariation(configPath("Epithelial", "cycle", "rate", 0), 0.0),
]

out = run(inputs, dvs...; n_replicates = 1, force_recompile = true)
sim_id = only(simulationIDs(out))
println("sim_id = ", sim_id)

outdir = mkpath(joinpath(@__DIR__, "figure3"))
snapshot_path(idx) = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output",
                               idx === :final ? "final.svg" : @sprintf("snapshot%08d.svg", idx))

MILESTONES = [125, 150, 175, 200]
for n in MILESTONES
    res = output_index_at_cell_count(sim_id, n)
    if res === nothing
        @warn "never reached $n cells — check the run"
        continue
    end
    idx, actual = res
    dst = joinpath(outdir, "n$(actual).svg")
    cp(snapshot_path(idx), dst; force = true)
    println("  -> ", dst, "  (target $n, actual $actual)")
end