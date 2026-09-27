using PhysiCellModelManager
using Distributions
include("GenerateReport.jl")   # pulls in ParameterOptimization.jl (summary + distance)

setNumberOfParallelSims(10)

# ────────────────────────────────────────────────────────────────────
# SANITY CHECK on an existing sim (run BEFORE resetDatabase wipes it)
# ────────────────────────────────────────────────────────────────────
println("\n=== Sanity check: trajectory + stability on sim 2 ===")

times, series = boundary_trajectory(2)
IC_t = getfield.(series, :IC)

println("frames: $(length(times))  (first=$(first(times)), last=$(last(times)))")
println("IC[start] = $(round(first(IC_t); digits=3))   IC[end] = $(round(last(IC_t); digits=3))")
println("IC time-course (every ~10th frame):")
for i in 1:10:length(times)
    println("  frame $(times[i]):  IC = $(round(IC_t[i]; digits=3))")
end

sm = stability_metrics(2; window = 20)
println("stability_metrics(2): terminal_drift = $(round(sm["terminal_drift"]; digits=5)), " *
        "settled IC = $(round(sm["IC"]; digits=3)), n_frames = $(Int(sm["n_frames"]))")
println("=== end sanity check ===\n")

println("=== Setting up ABC-SMC EQUILIBRIUM/STABILITY Calibration (Linear Mechanics) ===")
resetDatabase(; force_reset=true, force_continue=true)
println("Reset Database")

# ────────────────────────────────────────────────────────────────────
# Config
# ────────────────────────────────────────────────────────────────────
config_folder = "DuctDev_ParamOpt"
inputs = InputFolders(config_folder, config_folder; rulesets_collection=config_folder)

# Lock in the LINEAR mechanics (same regime as OptSMC.jl). Proliferation is
# driven by the ruleset in this config, so the duct still deforms — we are NOT
# turning it off; we want it to deform and then HOLD.
ref_model = createTrial(inputs,
    DiscreteVariation(configPath("max_time"), 11520), # 8 days

    # Strain Mechanics (ON, Linear)
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1),
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 1),

    # Restoring Mechanics (ON, Linear)
    DiscreteVariation(configPath("user_parameters", "is_restore_lin"), 1),

    # Lumenal Pressure (ON)
    DiscreteVariation(configPath("user_parameters", "is_lumenal_pressure"), 1),

    # Kernel Size (Static)
    DiscreteVariation(configPath("user_parameters", "membrane_force_smoothing_sigma"), 25.0)
)

# ────────────────────────────────────────────────────────────────────
# Summary statistic: trajectory-based stability at a fixed tail window.
# We wrap stability_summary_statistic in a closure that fixes `window`, exactly
# like OptSMC.jl wraps bm_summary_statistic to fix target_cell_count.
# ────────────────────────────────────────────────────────────────────
const STABILITY_WINDOW = 20   # ~last 20 boundary frames (~last 10–20% of the run)
summary_stability(monad_id) = stability_summary_statistic(monad_id; window = STABILITY_WINDOW)

observed_target = Dict(
    "terminal_drift" => 0.0,
    "IC"             => 1.41,   # ← REPLACE with the settled IC of your reference duct
)

# ────────────────────────────────────────────────────────────────────
# Search space (priors). Same mechanics knobs as OptSMC.jl.
# ────────────────────────────────────────────────────────────────────
parameters_to_tune = [
    DistributedVariation(configPath("user_parameters", "seg_lin"),  Uniform(0.01, 0.1)),
    DistributedVariation(configPath("user_parameters", "home_lin"), Uniform(0.0005, 0.005)),
    DistributedVariation(configPath("user_parameters", "lumenal_pressure_strength"), Uniform(0.005, 0.05)),
]

# ────────────────────────────────────────────────────────────────────
# Build & run calibration
# ────────────────────────────────────────────────────────────────────
println("Building Stability Calibration Problem...")
problem = CalibrationProblem(
    ref_model,
    parameters_to_tune,
    observed_target,
    summary_stability,       # ← trajectory-based stability summary (tail window)
    bm_distance;             # ← same scale-normalized distance; scales include terminal_drift
    n_replicates = 2
)

println("Running ABC-SMC (Stability run)...")
result = runABC(
    problem;
    population_size      = 50,
    max_nr_populations   = 6,
    minimum_epsilon      = 0.1,    # normalized-distance scale; tunable
    epsilon_quantile     = 0.5,
    min_acceptance_rate  = 0.02,   # stop if proposals stop getting accepted
    min_epsilon_decrease = 0.02,   # stop once ε plateaus
    description          = "Linear Mechanics EQUILIBRIUM/stability calibration",
)

# ────────────────────────────────────────────────────────────────────
# Extract & report
# ────────────────────────────────────────────────────────────────────
df, weights = posterior(result)
println("\n═══ Final Posterior Parameter Estimates (stable-and-deformed) ═══")
println(df)
println("\n═══ Convergence (ε / acceptance_rate / ess per generation) ═══")
println(ConvergenceSummary(result))

# posterior() returns only parameter columns — get sim IDs from the calibration's monads.
monad_ids = PhysiCellModelManager.ModelManager.calibrationMonadIDs(result.calibration)
calibration_sim_ids = sort!(unique!(reduce(vcat, simulationIDs.(Monad.(monad_ids)); init = Int[])))

println("\n═══ Generating Report for Best Fits ═══")
GenerateReport("Reports/Calibration_Results_Equilibrium";
    sim_ids = calibration_sim_ids,
    layout = :sequential,          # one snapshot per posterior parameter set, in order
    snapshot_selection = :final
)