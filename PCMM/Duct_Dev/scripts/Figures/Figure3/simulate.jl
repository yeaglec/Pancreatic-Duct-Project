# Figure 3 (simulations): equilibrium (no proliferation) vs the full default model (cancer cells proliferate).
# Edit: ROWS (name => parameter changes from the default model), TIMES (snapshot times, min), MAX_TIME.
# Output: Outputs/<row>/t<minutes>.svg — the same times in every row, so columns line up.
include(joinpath(@__DIR__, "..", "common.jl"))

MAX_TIME = 10080
TIMES = [0, 5040, 10080]   # day 0, 3.5, 7

ROWS = [
    ("equilibrium", [DiscreteVariation(configPath("CAF", "cycle", "rate", 0), 0.0)]),   # cancer-cell division off
    ("growth",      []),                                                                   # full default model
]

outdir_base = fresh_dir(joinpath(@__DIR__, "Outputs"))

for (i, (name, dvs)) in enumerate(ROWS)
    out = run(inputs, DiscreteVariation(configPath("max_time"), MAX_TIME), dvs...;
              n_replicates = 1, force_recompile = i == 1)   # both rows share one executable
    sim_id = only(simulationIDs(out))
    println(name, " -> sim_id ", sim_id)

    outdir = mkpath(joinpath(outdir_base, name))
    for t in TIMES
        save_frame(sim_id, frame_at(t), joinpath(outdir, "t$(t).svg"))
    end
end
