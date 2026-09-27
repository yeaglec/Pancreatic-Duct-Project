using PhysiCellModelManager
include("GenerateReport.jl")
setNumberOfParallelSims(10)
# Clear all previous simulations
# println("Clearing previous simulations database...")
# resetDatabase(; force_reset=true, force_continue=true)

############ set up ############

config_folder = custom_code_folder = "Duct_Dev_4_29" # this folder is located at Duct_Dev/data/inputs/configs
rulesets_collection_folder = "Duct_Dev_4_29" # this folder is located at Duct_Dev/data/inputs/rulesets_collections

# package them all together into a single object
inputs = InputFolders(config_folder, custom_code_folder;
                        rulesets_collection=rulesets_collection_folder,
)

# Reference variations
xml_path = configPath("max_time")
value = 10080
dv_max_time = DiscreteVariation(xml_path, value)

force_recompile = false
use_previous = true


############ run the sampling ############

# We will set the default simulations to have a lower max time.
# This will serve as a reference for the following simulations.

# EXPERIMENT 0 (Baseline)
# 6 threads
# println("Running Baseline...")
# n_replicates = 2
# out_baseline = run(inputs, dv_max_time; n_replicates=n_replicates, force_recompile=force_recompile)
# # makeMovie(out_baseline)

# # 1 thread recompile
# println("Running Baseline...")
# n_replicates = 2
# xml_path = configPath("max_time")
# value = 9000
# dv_max_time = DiscreteVariation(xml_path, value)
# out_baseline = run(inputs, dv_max_time; n_replicates=n_replicates, force_recompile=true)

# # 1 THREAD discrete var and recompile
# println("Running Baseline...")
# n_replicates = 2
# xml_path = configPath("max_time")
# value = 8000
# dv_max_time = DiscreteVariation(xml_path, value)
# out_baseline = run(inputs, dv_max_time; n_replicates=n_replicates, force_recompile=true)



# # EXPERIMENT 1
# println("Running Experiment 1...")

# xml_path_sigma = configPath("user_parameters", "membrane_force_smoothing_sigma")
# values_sigma = [ 25, 30, 34, 38, 42]                                 
# dv_sigma = DiscreteVariation(xml_path_sigma, values_sigma)

# out_exp1 = run(inputs, dv_sigma, dv_max_time; force_recompile=force_recompile)  #exp 2-7

# # EXPERIMENT 2
# println("Running Experiment 2...")
# xml_path_strain = configPath("user_parameters", "is_restore_lin")
# values_strain = [1]
# dv_strain = DiscreteVariation(xml_path_strain, values_strain)
# xml_path_sigma = configPath("user_parameters", "home_lin")
# values_sigma = [.0005,.00075, .001, .002, .0025, .0033] #exp 8-26
# dv_sigma = DiscreteVariation(xml_path_sigma, values_sigma)


# out_exp1 = run(inputs, dv_sigma, dv_strain,  dv_max_time; force_recompile=force_recompile)

# # EXPERIMENT 3
# println("Running Experiment 3...")
# xml_path_strain = configPath("user_parameters", "is_restore_lin")
# values_strain = [0]
# dv_strain = DiscreteVariation(xml_path_strain, values_strain)
# xml_path_sigma = configPath("user_parameters", "home_lin")
# values_sigma = [.002, .0025, .005]
# dv_sigma = DiscreteVariation(xml_path_sigma, values_sigma)
# xml_path_exp = configPath("user_parameters", "home_exp")
# values_exp= [.01, .02, .03]
# dv_exp = DiscreteVariation(xml_path_exp, values_exp)

# out_exp1 = run(inputs, dv_sigma, dv_strain, dv_exp, dv_max_time; force_recompile=force_recompile)

# GenerateReport("TEST2";
#     sim_ids = collect(7:15),
#     baseline_sim_id = 1,
#     snapshot_selection = :by_cell_count,
#     target_cell_count = 300,
# )

############ BM Evaluation Metrics ############
# These functions are defined in GenerateReport.jl.
# Each takes a simulation ID and an optional `timestep` keyword.
# If timestep is omitted, the final (largest) timestep is used.
#
#   Perimeter(sim_id)                — total edge length of BM polygon
#   Area(sim_id)                     — enclosed area (Shoelace formula)
#   AreaToPerimeterRatio(sim_id)     — A / P  (compactness)
#   InverseCircularity(sim_id)       — P²/(4πA), 1.0 = perfect circle
#   boundary_timesteps(sim_id)       — list available timestep indices
#
# Example: compare initial vs final membrane shape for simulation 1

# sim = 1
# ts = boundary_timesteps(sim)
# t0, tf = first(ts), last(ts)

# println("\n BM Metrics for Simulation $sim ")
# println("  Available timesteps: $(length(ts)) total")

# println("\n  t = $t0 (initial):")
# println("    Perimeter            = $(Perimeter(sim; timestep=t0))")
# println("    Area                 = $(Area(sim; timestep=t0))")
# println("    Area/Perimeter       = $(AreaToPerimeterRatio(sim; timestep=t0))")
# println("    Inverse Circularity  = $(InverseCircularity(sim; timestep=t0))")

# println("\n  t = $tf (final):")
# println("    Perimeter            = $(Perimeter(sim; timestep=tf))")
# println("    Area                 = $(Area(sim; timestep=tf))")
# println("    Area/Perimeter       = $(AreaToPerimeterRatio(sim; timestep=tf))")
# println("    Inverse Circularity  = $(InverseCircularity(sim; timestep=tf))")

# println("\n══ Done ══\n")

# EXPERIMENT 4
# println("Running Experiment 4...")
# xml_path_strain = configPath("user_parameters", "is_bending_stiffness")
# values_strain = [1]
# dv_strain = DiscreteVariation(xml_path_strain, values_strain)
# xml_path_sigma = configPath("user_parameters", "membrane_bending_constant")
# values_sigma = [.75, 1, 2, 3, 4, 5, 6, 8, 10]
# dv_sigma = DiscreteVariation(xml_path_sigma, values_sigma)

# out_exp1 = run(inputs, dv_sigma, dv_strain, dv_max_time; force_recompile=force_recompile)


# println("Running Experiment 5...")
# xml_path_strain = configPath("user_parameters", "is_add_membrane_nodes")
# values_strain = [1]
# dv_strain = DiscreteVariation(xml_path_strain, values_strain)
# xml_path_sigma = configPath("user_parameters", "max_edge_length")
# values_sigma = [2,4,8,10,15,20,25]
# dv_sigma = DiscreteVariation(xml_path_sigma, values_sigma)

# out_exp1 = run(inputs, dv_sigma, dv_strain, dv_max_time; force_recompile=force_recompile)


# println("Running Experiment 6...")
# xml_path_elast = configPath("user_parameters", "Segment_Elasticity")
# values_elast = [1]
# dv_elast = DiscreteVariation(xml_path_elast, values_elast)
# xml_path_strain = configPath("user_parameters", "is_strain_lin")
# values_strain = [0]
# dv_strain = DiscreteVariation(xml_path_strain, values_strain)
# xml_path_sigma = configPath("user_parameters", "seg_lin")
# values_sigma = [.01, .025, .05, .075, .1]
# dv_sigma = DiscreteVariation(xml_path_sigma, values_sigma)
# xml_path_exp = configPath("user_parameters", "seg_exp")
# values_exp= [.5, 1, 2, 4, 8]
# dv_exp = DiscreteVariation(xml_path_exp, values_exp)

# out_exp1 = run(inputs, dv_sigma, dv_strain, dv_elast, dv_exp, dv_max_time; force_recompile=force_recompile)


# println("Generating movies for Exp 1!")
# for id in 1:8
#     try
#         makeMovie(id)
#     catch e
#         println("Skipping movie for simulation \$id...")
#     end
# end

println("All experiments setup completely in GenerateData.jl! Ready to exit.")

# Okay, awesome all my simulations are finished. Now I want you to go through each simulation (2-71) and take the snapshot00001200.svg file from each of their output folders (PCMM/Duct_Dev/data/outputs/simulations/"the sim number"/output/snapshot00001200.svg) and move it into the new folder I just created called Model Comparisons. Here is the tricky part, I want you to rename the snapshot00001200.svg according to the parameter being changed in the generatedata script. For example, simulation 2 set "membrane_force_smoothing_sigma" to 8 and simulation 3 set sigma to 16. So you would call the snapshot0001200.svg file for simulation 2 sigma_8 for example. Same for experiment 2. For example simulation 8 call those snapshot linrestore_.0005)


############ ABC-SMC Calibration ############
# ParameterOptimization.jl is included via GenerateReport.jl
using Distributions

resetDatabase(; force_reset=true, force_continue=true)


# Fix non-calibrated parameters via a reference monad (n_replicates=0  no sims run)

println("Beginning Optimization:")
ref = createTrial(inputs,
    DiscreteVariation(configPath("user_parameter","is_gaussian_smoothing"), 1),
    DiscreteVariation(configPath("user_parameter","Segment_Elasticity"), 0),
    DiscreteVariation(configPath("max_time"), 10080);
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

# Build calibration problem
problem = CalibrationProblem(
    ref,                     # Monad 
    parameters,              # DistributedVariations with priors
    observed,                # target summary statistics
    bm_summary_statistic,    # (monad_id → Dict) 
    bm_distance;             # (sim, obs → Float64) 
    n_replicates = 2,
)

# resumeABC(Calibration(id))
# Run ABC-SMC — start small to verify the pipeline works
result = runABC(
    problem;
    population_size    = 10,    # small for testing; increase to 50-200 for real runs
    max_nr_populations = 3,     # small for testing; increase to 10-15 for real runs
    minimum_epsilon    = 0.05,
    description        = "BM morphology calibration (test run)",
)

# Extract the posterior distribution of best-fit parameters
df, weights = posterior(result)
println("\n═══ Posterior Parameter Estimates ═══")
println(df)

############ Report Generation ############
# Visualize the ABC-SMC calibration results
# Extracts sim IDs from the posterior and generates a comparison report
calibration_sim_ids = sort(collect(df.sim_id))
println("\n═══ Generating Report for Calibration Results ═══")
GenerateReport("Calibration_Report";
    sim_ids = calibration_sim_ids,
    snapshot_selection = :final,
)