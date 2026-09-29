# Figure 1 (simulation part): one cancer cell in the default model, to show extrusion after its first division.
# Edit: dvs (run settings) and the exploration window (k range, in 6-min frames around the division).
# Output: Outputs/explore/k<offset>.svg — every frame around the first division; compose.jl picks from these.
include(joinpath(@__DIR__, "..", "common.jl"))

outdir = fresh_dir(joinpath(@__DIR__, "Outputs"))

# One cancer cell in the ring; all mechanics at config defaults (only cancer cells divide by default)
dvs = [
    DiscreteVariation(configPath("max_time"), 2880),
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

# EXPLORATION: every frame from 1 before the division to 3 h after (6 min saves) — pick the extrusion frames from these
explore = mkpath(joinpath(outdir, "explore"))
for k in -1:30
    src = snapshot_path(sim_id, div_idx + k)
    isfile(src) || break
    cp(src, joinpath(explore, @sprintf("k%+03d.svg", k)); force = true)
end
println("  -> ", explore)
