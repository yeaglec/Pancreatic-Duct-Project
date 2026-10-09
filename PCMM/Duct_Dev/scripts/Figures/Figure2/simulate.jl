# Figure 2 (simulations): a healthy duct at equilibrium (no cancer cells, epithelial proliferation on), and the
# same duct after each perturbation, shown separately.
# Edit: BASELINE (shared settings), ROWS (name => the perturbation it adds), MAX_TIME.
# Output: Outputs/<row>/h<hours>.svg (one frame per hour) and Outputs/<row>/counts.csv (time in min, total cells)
include(joinpath(@__DIR__, "..", "common.jl"))

MAX_TIME = 10080   # 7 days
BASELINE = [
    DiscreteVariation(configPath("max_time"), MAX_TIME),
    DiscreteVariation(configPath("user_parameters", "num_cancer_cells"), 0),
    DiscreteVariation(configPath("Epithelial", "cycle", "rate", 0), 0.00144),   # same as the cancer cells; XML default is 0 (no EP division)
    DiscreteVariation(configPath("user_parameters", "neighbor_contact_inhibition"), 1),   # divide only with a free side (see Figures/README)
]

ROWS = [
    ("control", []),
    ("gap", [DiscreteVariation(configPath("user_parameters", "num_removed_cells"), 5)]),
    ("displaced", [DiscreteVariation(configPath("user_parameters", "num_displaced_cells"), 5),
                   DiscreteVariation(configPath("user_parameters", "displaced_distance"), 25.0),              # ~27 µm surface gap
                   DiscreteVariation(configPath("user_parameters", "membrane_interaction_length"), 30.0)]),   # reach them (default 10)
]

outdir_base = fresh_dir(joinpath(@__DIR__, "Outputs"))

for (i, (name, dvs)) in enumerate(ROWS)
    out = run(inputs, BASELINE..., dvs...; n_replicates = 1, force_recompile = i == 1)   # all rows share one executable
    sim_id = only(simulationIDs(out))
    println(name, " -> sim_id ", sim_id)

    outdir = mkpath(joinpath(outdir_base, name))
    for h in 0:(MAX_TIME ÷ 60)
        save_frame(sim_id, frame_at(60h), joinpath(outdir, "h$(h).svg"))
    end

    # Total cell count at every saved frame, for the recovery trace
    spts = PhysiCellModelManager.SimulationPopulationTimeSeries(sim_id)
    open(joinpath(outdir, "counts.csv"), "w") do io
        println(io, "time,cells")
        for (k, t) in enumerate(spts.time)
            println(io, t, ",", sum(c[k] for c in values(spts.cell_count); init = 0))
        end
    end
end
