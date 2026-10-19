# Contact-inhibition test for Figure 2: does a healthy duct stay at equilibrium, and does a gap refill?
# Usage (from scripts/):  PCMM_NUM_PARALLEL_SIMS=6 julia --project=. Model_Evaluation/Contact_Inhibition_Tests/sweep.jl
# Edit: SETTINGS (Epithelial pressure rule: half_max, hill_power), BASELINE, ROWS.
# Output: Outputs/<setting>_<row>/h<hours>.svg (every 6 h) and counts.csv (time in min, total cells)
include(joinpath(@__DIR__, "..", "..", "Figures", "common.jl"))

MAX_TIME = 4320   # 3 days
BASELINE = [
    DiscreteVariation(configPath("max_time"), MAX_TIME),
    DiscreteVariation(configPath("user_parameters", "num_cancer_cells"), 0),
    DiscreteVariation(configPath("Epithelial", "cycle", "rate", 0), 0.00144),   # same as the cancer cells
]

SETTINGS = [(0.5, 4.0), (1.0, 10.0)]                 # (half_max, hill_power); (0.5, 4) is the current rule
ROWS = [("control", 0), ("gap", 5)]                  # (name, num_removed_cells)

rule(p) = rulePath("Epithelial", "cycle entry", "decreasing_signals", "signal:name:pressure", p)
outdir_base = fresh_dir(joinpath(@__DIR__, "Outputs"))

for (i, (half_max, hill)) in enumerate(SETTINGS)
    dvs = [BASELINE...,
           DiscreteVariation(rule("half_max"), half_max),
           DiscreteVariation(rule("hill_power"), hill),
           DiscreteVariation(configPath("user_parameters", "num_removed_cells"), last.(ROWS))]   # both rows in parallel
    out = run(inputs, dvs...; n_replicates = 1, force_recompile = i == 1)   # NOTE: PCMM reuses matching sims — reset the database after changing model code
    for sim_id in simulationIDs(out)
        # tell the rows apart by their starting cell count (the gap row starts with fewer cells)
        n0 = last(output_index_at_cell_count(sim_id, 0))
        name = n0 < 125 ? "gap" : "control"   # 125 = number_total_cells
        tag = "hm$(half_max)_hp$(Int(hill))_$(name)"
        println(tag, " -> sim_id ", sim_id)

        outdir = mkpath(joinpath(outdir_base, tag))
        for h in 0:6:(MAX_TIME ÷ 60)
            save_frame(sim_id, frame_at(60h), joinpath(outdir, "h$(h).svg"))
        end
        spts = PhysiCellModelManager.SimulationPopulationTimeSeries(sim_id)
        open(joinpath(outdir, "counts.csv"), "w") do io
            println(io, "time,cells")
            for (k, t) in enumerate(spts.time)
                println(io, t, ",", sum(c[k] for c in values(spts.cell_count); init = 0))
            end
        end
    end
end
