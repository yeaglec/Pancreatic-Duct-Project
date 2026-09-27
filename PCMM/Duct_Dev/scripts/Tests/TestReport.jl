#=
TestReport.jl — fast end-to-end check of GenerateReport.jl
=============================================================================
Mirrors ParamOpt.jl's Experiment 1 (runs 1A + 1B) but with max_time = 12 so each
simulation finishes ~instantly. Verifies that:
  • each run() becomes its own correctly-labeled section,
  • layout autodetection picks 1D (run 1A) vs. sequential-scattered-2D (run 1B),
  • the :sequential and :grid overrides work,
  • output is organized as <dir>/index.html + <dir>/snapshots/.

Non-destructive: it does NOT reset the database — it only appends the small test
runs and reports on their own simulation IDs. Run it from the PCMM project root:

    using PhysiCellModelManager
    include("Duct_Dev/scripts/TestReport.jl")
=#

using PhysiCellModelManager
using Distributions
include("GenerateReport.jl")

# Initialize + resolve report paths from this script's location (cwd-independent).
project_dir = dirname(@__DIR__)                       # .../PCMM/Duct_Dev
report_base = joinpath(dirname(project_dir), "Reports", "TestReports")  # .../PCMM/Reports/TestReports
initializeModelManager(project_dir)

setNumberOfParallelSims(10)

############ Setup ############
config_folder = custom_code_folder = rulesets_collection_folder = "DuctDev_ParamOpt"
inputs = InputFolders(config_folder, custom_code_folder; rulesets_collection=rulesets_collection_folder)

dv_max_time = DiscreteVariation(configPath("max_time"), 12)   # 12 min → instant sims
force_recompile = false

############ Run 1A: Linear strain (1 swept param → expect 1D) ############
println("Running 1A: Linear strain sweep...")
ref_strain_lin = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1),
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 1)
)
out_strain_lin = run(LHSVariation(5), ref_strain_lin,
    DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.0005, 0.01));
    n_replicates = 1, force_recompile = force_recompile
)

############ Run 1B: Exponential strain (2 swept params via LHS → expect sequential) ############
println("Running 1B: Exponential strain sweep...")
ref_strain_exp = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1),
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 0)
)
out_strain_exp = run(LHSVariation(6), ref_strain_exp,
    DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.01, 0.1)),
    DistributedVariation(configPath("user_parameters", "seg_exp"), Uniform(0.01, 0.05));
    n_replicates = 1, force_recompile = force_recompile
)

############ Generate test reports ############
println("\n=== TEST 1: sections + :auto  → two sections (Linear=1D, Exponential=sequential) ===")
GenerateReport(joinpath(report_base, "Test_Sections");
    sections = ["Linear Strain"      => out_strain_lin,
                "Exponential Strain" => out_strain_exp],
    snapshot_selection = :final
)

println("\n=== TEST 2: flat sim_ids + layout=:sequential  → one sequential section (OptSMC-style) ===")
GenerateReport(joinpath(report_base, "Test_Sequential");
    sim_ids = vcat(simulationIDs(out_strain_lin), simulationIDs(out_strain_exp)),
    layout = :sequential,
    snapshot_selection = :final
)

println("\n=== TEST 3: layout=:grid  → forced (sparse) 2D tableau on the 2 swept params ===")
GenerateReport(joinpath(report_base, "Test_Grid");
    sections = ["Exponential Strain" => out_strain_exp],
    layout = :grid,
    snapshot_selection = :final
)

println("\nAll test reports generated under: $(report_base)")
