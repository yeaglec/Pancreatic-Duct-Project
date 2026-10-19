# Runs every experiment for the baseline and one candidate (1 simulation each) and records which simulation is which.
# Usage (from scripts/):  PCMM_NUM_PARALLEL_SIMS=4 julia --project=. Model_Testing/run_comparison.jl <candidate>
#   <candidate> = a file name in Candidates/ without .jl, e.g. adhesion_0.1
#   MAX_TIME=60 in the shell shortens every experiment (a quick check that everything runs)
# Output: Outputs/<candidate>/simulations.csv (experiment, version, inputs, sim_id) and a copy of the candidate file
using PhysiCellModelManager

initializeModelManager(joinpath(@__DIR__, "..", ".."))   # -> Duct_Dev
include(joinpath(@__DIR__, "experiments.jl"))

BASELINE = "Duct_Baseline_V0.0.1"
SAVE_INTERVAL = 30   # min between SVG / full-data saves (does not change the simulation itself)

isempty(ARGS) && error("name a candidate, e.g.  julia --project=. Model_Testing/run_comparison.jl adhesion_0.1")
candidate_file = joinpath(@__DIR__, "Candidates", ARGS[1] * ".jl")
candidate = include(candidate_file)
experiments = load_experiments()

folders(name) = InputFolders(name, name; rulesets_collection = name)
run_length(e) = haskey(ENV, "MAX_TIME") ? parse(Int, ENV["MAX_TIME"]) : e.max_time
settings(e) = [DiscreteVariation(configPath("max_time"), run_length(e)),
               DiscreteVariation(configPath("SVG_save_interval"), SAVE_INTERVAL),
               DiscreteVariation(configPath("full_data_interval"), SAVE_INTERVAL)]

# One simulation per experiment and version. PCMM reuses a stored simulation when everything matches,
# so the baseline runs only once per experiment across all comparisons.
baseline_runs = [createTrial(folders(BASELINE), settings(e)..., e.changes...; n_replicates = 1) for e in experiments]
candidate_runs = [createTrial(folders(candidate.inputs), settings(e)..., e.changes..., candidate.changes...; n_replicates = 1)
                  for e in experiments]
run(Trial([Sampling(baseline_runs), Sampling(candidate_runs)]))   # all in parallel

outdir = mkpath(joinpath(@__DIR__, "Outputs", candidate.name))
cp(candidate_file, joinpath(outdir, "candidate.jl"); force = true)
open(joinpath(outdir, "simulations.csv"), "w") do io
    println(io, "experiment,version,inputs,sim_id")
    for (e, b, c) in zip(experiments, baseline_runs, candidate_runs)
        println(io, e.name, ",baseline,", BASELINE, ",", only(simulationIDs(b)))
        println(io, e.name, ",candidate,", candidate.inputs, ",", only(simulationIDs(c)))
    end
end
println(read(joinpath(outdir, "simulations.csv"), String))
