# Figure 4 (simulations): sweep each active membrane mechanic around its default, compared at one growth milestone.
# Edit: SWEEPS (name, XML path, values — include the default as a reference), TARGET_CELL_COUNT, MAX_TIME.
# Output: Outputs/<name>/v<value>.svg — the first frame reaching TARGET_CELL_COUNT (final frame if never reached).
include(joinpath(@__DIR__, "..", "common.jl"))

MAX_TIME = 10080
TARGET_CELL_COUNT = 175   # compare all variants at the same growth milestone

# Defaults (Duct_Dev_Figs): seg_lin 0.1, membrane_bending_constant 4, lumenal_pressure_strength 0.01
SWEEPS = [
    ("strain",   configPath("user_parameters", "seg_lin"),                   [0.025, 0.05, 0.1, 0.2, 0.4]),
    ("bending",  configPath("user_parameters", "membrane_bending_constant"), [0.5, 1.0, 2.0, 4.0, 8.0]),
    ("pressure", configPath("user_parameters", "lumenal_pressure_strength"), [0.0025, 0.005, 0.01, 0.02, 0.04]),
]

outdir_base = fresh_dir(joinpath(@__DIR__, "Outputs"))

first_run = true
for (name, path, values) in SWEEPS
    outdir = mkpath(joinpath(outdir_base, name))
    for v in values
        # the default value repeats across sweeps; PCMM reuses that simulation instead of rerunning it
        out = run(inputs, DiscreteVariation(configPath("max_time"), MAX_TIME), DiscreteVariation(path, v);
                  n_replicates = 1, force_recompile = first_run)
        global first_run = false
        sim_id = only(simulationIDs(out))
        println(name, " = ", v, " -> sim_id ", sim_id)
        save_snapshot(sim_id, TARGET_CELL_COUNT, joinpath(outdir, "v$(v).svg"))
    end
end
