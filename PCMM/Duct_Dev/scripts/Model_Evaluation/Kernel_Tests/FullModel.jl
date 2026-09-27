#=

==============Important Note==============
This script is intedned to test the gaussian kernel and membrane strain mechanics on the FULL model.
We sweep throught the gaussian kernel and membrane strain parameters for testing 



================================================================================
PCMM SYNTAX REFERENCE & EXAMPLES
================================================================================

PARAMETER VARIATIONS
   - DiscreteVariation (Fixed parameter value):
     DiscreteVariation(configPath("group", "param"), value)
     Example: DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1)

   - DistributedVariation (Range for parameter sweeps):
     DistributedVariation(configPath("group", "param"), Distribution)
     Example: DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.01, 0.2))

LHS PARAMETER SWEEPS
   Runs simulations sampled across the parameter space.
   
   Syntax:
     run(inputs, variations...; sampler=LatinHypercube(N), n_replicates=R)
   Example:
     out = run(inputs, dv_max_time, 
         DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.01, 0.2));
         sampler = LatinHypercube(25), n_replicates = 1
     )

COMPARATIVE REPORTS
   Generates interactive HTML dashboards comparing simulation metrics.
   
   Syntax:
     GenerateReport("Report_Name"; sim_ids=Vector{Int}, kwargs...)
     * snapshot_selection options: :final, :by_cell_count, :by_index
   Example:
     GenerateReport("Report_1_Strain_Force";
         sim_ids = simulationIDs(out),
         snapshot_selection = :by_cell_count,
         target_cell_count = 200
     )

================================================================================
=#

using PhysiCellModelManager
using Distributions # Required for Uniform() sampling
include(joinpath(@__DIR__, "..", "GenerateReport.jl"))

setNumberOfParallelSims(10)

# Clear all previous simulations to ensure a clean slate for the reports
println("Clearing previous simulations database...")
resetDatabase(; force_reset=true, force_continue=true)

############ Setup ############

config_folder = custom_code_folder = rulesets_collection_folder = "DuctDev_ParamOpt" 
inputs = InputFolders(config_folder, custom_code_folder;
                      rulesets_collection=rulesets_collection_folder)

# Set max time to 5 days 
dv_max_time = DiscreteVariation(configPath("max_time"), 17280)
force_recompile = false

############ EXPERIMENT 1: Strain Only Sweep ############
println("Running Experiment 1: Strain Force (Linear & Exponential)...")

# 1A: Linear Strain Sweep
ref_strain_lin = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1), # Turn ON Strain
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 1),      # Use Linear
    DiscreteVariation(configPath("user_parameters", "is_gaussian_smoothing"), 0)      # Use Linear
)
out_strain_lin = run(LHSVariation(10), ref_strain_lin,
    DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.05, 0.75));
    n_replicates = 1,
    force_recompile = force_recompile
)

############ EXPERIMENT 2: Kernel only Sweep ############

# 1A: Linear Strain Sweep
ref_kernel = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 0), # Turn off Strain
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 0),     
    DiscreteVariation(configPath("user_parameters", "is_gaussian_smoothing"), 1)      # Mkae sure to turn on kernel


)
out_kernel = run(LHSVariation(10), ref_kernel,
    DistributedVariation(configPath("user_parameters", "membrane_force_smoothing_sigma"), Uniform(3, 20));
    n_replicates = 1,
    force_recompile = force_recompile
)

############ EXPERIMENT 3: Kernel and Strain Sweep ############
println("Running Experiment 3: Kernel and Strain...")

# 3A: Kernel and Strain Sweep
ref_kernel_strain = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1), # Turn on Strain
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 1),      # Use Linear
    DiscreteVariation(configPath("user_parameters", "is_gaussian_smoothing"), 1)      # Mkae sure to turn on kernel


)
out_kernel_strain = run(LHSVariation(20), ref_kernel_strain,
    DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.05, 0.75)),
    DistributedVariation(configPath("user_parameters", "membrane_force_smoothing_sigma"), Uniform(3, 20));
    n_replicates = 1,
    force_recompile = force_recompile
)

############ EXPERIMENT 4: No Kernel, Strain, and Bending  ############
println("Running Experiment 4: No Kernel, Strain, and Bending...")

# 4A: Linear Restoring Force Sweep
ref_strain_bending = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1), # Turn on Strain
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 1),      # Use Linear
    DiscreteVariation(configPath("user_parameters", "is_gaussian_smoothing"), 0),     # Mkae sure to turn on kernel
    DiscreteVariation(configPath("user_parameters", "is_bending_stiffness"), 1)      # Mkae sure to turn on kernel
)
out_strain_bending = run(LHSVariation(20), ref_strain_bending,
    DistributedVariation(configPath("user_parameters", "membrane_bending_constant"), Uniform(.1, 10)),
    DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.05, 0.75));
    n_replicates = 1,
    force_recompile = force_recompile
)

############ EXPERIMENT 5: No Kernel, Strain, and Bending  ############
println("Running Experiment 5: Kernel, Strain, and Bending...")

# 4A: Linear Restoring Force Sweep
ref_strain_bending_kernel = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1), # Turn on Strain
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 1),     # Use Linear
    DiscreteVariation(configPath("user_parameters", "is_gaussian_smoothing"), 1),      # Mkae sure to turn on kernel
    DiscreteVariation(configPath("user_parameters", "membrane_force_smoothing_sigma"), 3),     # Mkae sure to turn on kernel
    DiscreteVariation(configPath("user_parameters", "is_bending_stiffness"), 1)      # Mkae sure to turn on kernel

)
out_strain_bending_kernel = run(LHSVariation(20), ref_strain_bending_kernel,
    DistributedVariation(configPath("user_parameters", "membrane_bending_constant"), Uniform(.1, 10)),
    DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.05, 0.75));
    n_replicates = 1,
    force_recompile = force_recompile
)


############ Generate Reports ############
println("Simulations complete. Generating specific reports...")

GenerateReport("Reports/Kernel_Tests/Report_1_Kernel_Strain";
    sections = ["Kernel Only" => out_kernel,
                "Strain Only" => out_strain_lin,
                "Kernel and Strain" => out_kernel_strain],
    snapshot_selection = :by_cell_count,
    target_cell_count = 200
)

GenerateReport("Reports/Kernel_Tests/Report_2_Bending_Force";
    sections = [
                "No Kernel, Strain, and Bending" => out_strain_bending,
                "Kernel, Strain, and Bending" => out_strain_bending_kernel],
    snapshot_selection = :by_cell_count,
    target_cell_count = 200
)


println("All parameter sweeps and reports successfully generated!")