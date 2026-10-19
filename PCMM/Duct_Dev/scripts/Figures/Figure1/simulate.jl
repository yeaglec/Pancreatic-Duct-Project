# Figure 1 (simulation part): one cancer cell in the default model, from its first division to a second layer.
# Edit: dvs (run settings) and the exploration frames (fine around the division, hourly over the whole run).
# Output: Outputs/explore/k<offset>.svg — offset in 6-min frames from the first division; compose.jl picks from these.
include(joinpath(@__DIR__, "..", "common.jl"))

outdir = fresh_dir(joinpath(@__DIR__, "Outputs"))

# One cancer cell in the ring; all mechanics at config defaults (only cancer cells divide by default)
dvs = [
    DiscreteVariation(configPath("max_time"), 5760),   # 4 days
    DiscreteVariation(configPath("user_parameters", "num_cancer_cells"), 1),
]

out = run(inputs, dvs...; n_replicates = 1, force_recompile = true)
sim_id = only(simulationIDs(out))
println("sim_id = ", sim_id)

# First division = first frame with one more cell than at the start
n0 = last(output_index_at_cell_count(sim_id, 0))   # frame 0 always qualifies -> initial cell count
div_idx = snapshot_index(sim_id, n0 + 1)
div_idx === :final && error("the cancer cell never divided — increase max_time")
println("first division at snapshot ", div_idx)

# EXPLORATION: every frame for 3 h after the division (6 min saves) + one frame per hour over the whole run
explore = mkpath(joinpath(outdir, "explore"))
n_frames = count(f -> startswith(f, "snapshot"), readdir(dirname(snapshot_path(sim_id, 0))))
fine = -1:30
hourly = (-div_idx:frame_at(60):(n_frames - 1 - div_idx))
for k in sort(union(fine, hourly))
    src = snapshot_path(sim_id, div_idx + k)
    isfile(src) && cp(src, joinpath(explore, @sprintf("k%+04d.svg", k)); force = true)
end
println("  -> ", explore)
