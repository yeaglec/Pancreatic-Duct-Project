# Figure 2 (simulations): the default model in each initial membrane geometry.
# Edit: SHAPES (membrane_shape value => name), TIMES (snapshot times, min), MAX_TIME, and BASELINE.
# Output: Outputs/<shape>/t<minutes>.svg
include(joinpath(@__DIR__, "..", "common.jl"))

MAX_TIME = 10080
BASELINE = [
    DiscreteVariation(configPath("max_time"), MAX_TIME),
    DiscreteVariation(configPath("user_parameters", "membrane_circle_radius"), 275.0),   # circle only: match the default shape's size (R_eff ~ 274 µm; XML default is 300)
]

SHAPES = [(0, "default"), (1, "star"), (2, "circle")]   # (membrane_shape value, name)
TIMES = [0, 3360, 6720, 10080]                            # 4 equally spaced time points (min)

outdir_base = fresh_dir(joinpath(@__DIR__, "Outputs"))

for (i, (value, name)) in enumerate(SHAPES)
    dvs = [BASELINE..., DiscreteVariation(configPath("user_parameters", "membrane_shape"), value)]
    out = run(inputs, dvs...; n_replicates = 1, force_recompile = i == 1)   # all shapes share one executable
    sim_id = only(simulationIDs(out))
    println(name, " -> sim_id ", sim_id)

    outdir = mkpath(joinpath(outdir_base, name))
    for t in TIMES
        save_frame(sim_id, frame_at(t), joinpath(outdir, "t$(t).svg"))
    end
end
