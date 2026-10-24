# Reference run of the current baseline: default settings, 7 days, one simulation.
# Usage (from scripts/):  julia --project=. Model_Testing/run_baseline.jl
# Output: Outputs/baseline_<version>_final.svg (the full outputs stay in data/outputs/simulations/<id>)
using PhysiCellModelManager

initializeModelManager(joinpath(@__DIR__, "..", ".."))   # -> Duct_Dev

BASELINE = "Duct_Baseline_V0.0.1"
MAX_TIME = 10080   # 7 days
setNumberOfParallelSims(4)

inputs = InputFolders(BASELINE, BASELINE; rulesets_collection = BASELINE)
out = run(inputs, DiscreteVariation(configPath("max_time"), MAX_TIME); n_replicates = 1, force_recompile = true)
sim_id = only(simulationIDs(out))

sim_dir = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output")
outdir = mkpath(joinpath(@__DIR__, "Outputs"))
cp(joinpath(sim_dir, "final.svg"), joinpath(outdir, "baseline_$(BASELINE)_final.svg"); force = true)
println(BASELINE, " -> sim_id ", sim_id)
