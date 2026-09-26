using PhysiCellModelManager
include("GenerateReport.jl")
setNumberOfParallelSims(10)

############ set up ############

config_folder = custom_code_folder = "Duct_Dev_4_29" # this folder is located at Duct_Dev/data/inputs/configs
rulesets_collection_folder = "Duct_Dev_4_29" # this folder is located at Duct_Dev/data/inputs/rulesets_collections

# package them all together into a single object
inputs = InputFolders(config_folder, custom_code_folder;
                        rulesets_collection=rulesets_collection_folder,
)

############ ABC-SMC Calibration (Fast Test) ############
# ParameterOptimization.jl is included via GenerateReport.jl
using Distributions

resetDatabase(; force_reset=true, force_continue=true)

# Fix non-calibrated parameters via a reference monad (n_replicates=0  no sims run)
# FAST TEST: max_time is reduced to 24 (minutes of simulation time) to run super fast.
println("Beginning Optimization:")
ref = createTrial(inputs,
    DiscreteVariation(configPath("user_parameter","is_gaussian_smoothing"), 1),
    DiscreteVariation(configPath("user_parameter","Segment_Elasticity"), 0),
    DiscreteVariation(configPath("max_time"), 1440); # Reduced from 10080 to 24 for fast test
    n_replicates = 0
)

println("DiscreteVariations!!!")
# Parameters to infer, with uniform priors
parameters = [
    DistributedVariation(configPath("user_parameters","seg_lin"), Uniform(0.01, 0.1)),
    DistributedVariation(configPath("user_parameters","seg_exp"), Uniform(0.5, 8.0)),
    DistributedVariation(configPath("user_parameters","home_lin"), Uniform(0.0005, 0.01)),
    DistributedVariation(configPath("user_parameters","home_exp"), Uniform(0.005, 0.05)),
]

# Observed (target) BM morphology — replace with ground truth stats when available
observed = Dict(
    "IC"                => 1.15,   
    "area_frac_change"  => -0.15,   
    "max_displacement"  => 40.0,   
)

# Step 1: First step determine parameters and get rough disterbutions (up or down 20%, maybe 50%)
# Step 2: soblo sequences are deterministic unlike LHS 
# Step 3: calibration
# 

# Build calibration problem
# FAST TEST: n_replicates is reduced to 1 for speed.
problem = CalibrationProblem(
    ref,                     # Monad 
    parameters,              # DistributedVariations with priors
    observed,                # target summary statistics
    bm_summary_statistic,    # (monad_id → Dict) 
    bm_distance;             # (sim, obs → Float64) 
    n_replicates = 1,        # Reduced from 2 to 1 for fast test
)

# Run ABC-SMC
# FAST TEST: population_size is reduced to 3 and max_nr_populations to 2.
result = runABC(
    problem;
    population_size    = 3,     # Reduced from 10 to 3 for fast test
    max_nr_populations = 2,     # Reduced from 3 to 2 for fast test
    minimum_epsilon    = 0.05,
    description        = "BM morphology calibration (fast test run)",
)

# Extract the posterior distribution of best-fit parameters
df, weights = posterior(result)
println("\n═══ Posterior Parameter Estimates ═══")
println(df)

############ Report Generation ############
# Visualize the ABC-SMC calibration results (fast test)
calibration_sim_ids = Int[]
for mid in df.monad_id
    append!(calibration_sim_ids, simulationIDs(Monad(mid)))
end
sort!(unique!(calibration_sim_ids))

println("\n═══ Generating Report for Calibration Results ═══")
GenerateReport("Reports/Calibration_Report_Fast";
    sim_ids = calibration_sim_ids,
    snapshot_selection = :final,
)
