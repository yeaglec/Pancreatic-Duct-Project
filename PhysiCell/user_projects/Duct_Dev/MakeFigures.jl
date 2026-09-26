#=
================================================================================
MakeFigures.jl — simulation snapshots for Results Claims 1–3
================================================================================

Deliberately minimal. This script does exactly one thing: run the handful of
simulations the first three claims need and copy the relevant PhysiCell SVG
snapshots into ./figures/ under readable names. No plotting libraries, no
compositing — you lay the panels out in LaTeX with subfigure.

    julia MakeFigures.jl

Outputs (figures/):
  claim1_rep{1,2,3}_{initial,final}.svg   equilibrium, no proliferation (Fig 1C)\n  claim1_equilibrium_rep{1,2,3}.csv       drift + IC per frame      (Fig 1D)
  claim3_n{125,150,175,200}.svg           proliferation -> inward deformation
  claim3_deformation.csv                  descriptor vs. cell count (Fig 3D)
================================================================================
=#

using PhysiCellModelManager
include("GenerateReport.jl")   # -> ParameterOptimization.jl helpers
                               #    (output_index_at_cell_count,
                               #     boundary_timestep_at_cell_count,
                               #     InverseCircularity, MaxNodeDisplacement)
                               #    and snapshot_svg_path

setNumberOfParallelSims(10)

# ───────────────────────────── config ─────────────────────────────
const FIG_DIR = "figures"
mkpath(FIG_DIR)

config_folder = custom_code_folder = rulesets_collection_folder = "DuctDev_ParamOpt"
inputs = InputFolders(config_folder, custom_code_folder;
                      rulesets_collection = rulesets_collection_folder)

dv_max_time     = DiscreteVariation(configPath("max_time"), 7200)   # 5 days
force_recompile = false

# The "full mechanics" configuration the figures are meant to show.
# Re-check these against your latest calibration before treating figures as final.
FULL_MECHANICS = [
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1),
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"),      1),
    DiscreteVariation(configPath("user_parameters", "is_restore_lin"),     1),
    DiscreteVariation(configPath("user_parameters", "is_lumenal_pressure"),1),
]
NO_PROLIFERATION = DiscreteVariation(configPath("user_parameters", "proliferation_exit_rate"), 0.0)

# ───────────────────────────── helpers ────────────────────────────

"Copy one PhysiCell snapshot SVG into figures/ under a readable name."
function save_snapshot(sim_id::Int, index, name::AbstractString)
    src = snapshot_svg_path(sim_id, index)
    isfile(src) || error("no snapshot for sim $sim_id at index $index: $src")
    dst = joinpath(FIG_DIR, string(name, ".svg"))
    cp(src, dst; force = true)
    println("  -> ", dst)
    return dst
end

"Frame index (and the cell count actually reached) nearest a cell-count milestone."
function frame_at(sim_id::Int, n::Int)
    res = output_index_at_cell_count(sim_id, n)
    res === nothing && error("sim $sim_id never reached $n cells — lengthen dv_max_time")
    return res   # (index, actual_count)
end

# ══════════════ CLAIM 1 — monolayer is stable at equilibrium ══════════════
# Same mechanics, proliferation off, three replicates: nothing should move.
# Panels C (snapshots) and D (drift traces) both come from these runs.
println("Claim 1: equilibrium (no proliferation)")

out1  = run(inputs, dv_max_time, FULL_MECHANICS..., NO_PROLIFERATION;
            n_replicates = 3, force_recompile = force_recompile)
sims1 = simulationIDs(out1)

# --- Panel C: first/last snapshot of each replicate -----------------------
for (k, sid) in enumerate(sims1)
    save_snapshot(sid, 0,      "claim1_rep$(k)_initial")
    save_snapshot(sid, :final, "claim1_rep$(k)_final")
end

# --- Panel D: does anything drift? ---------------------------------------
# Every boundary frame of every replicate, two numbers:
#   max_drift  max |node - its initial position|  (um)   -> BM should not move
#   IC         perimeter^2 / (4*pi*area)                 -> should sit at ~1
# Flat traces here are the noise floor that later deformation is measured against.

"Simulation time (min) per boundary frame; falls back to frame index."
function frame_times(sim_id::Int, ts::Vector{Int})
    try
        t = PhysiCellModelManager.SimulationPopulationTimeSeries(sim_id).time
        return [i + 1 <= length(t) ? t[i+1] : NaN for i in ts]   # ts is 0-based
    catch
        return Float64.(ts)
    end
end

# One file per replicate keeps the LaTeX side trivial (no row filtering).
for (k, sid) in enumerate(sims1)
    ts = boundary_timesteps(sid)
    tm = frame_times(sid, ts)
    path = joinpath(FIG_DIR, "claim1_equilibrium_rep$(k).csv")
    open(path, "w") do io
        println(io, "frame,time_h,max_drift,inverse_circularity")
        for (j, t) in enumerate(ts)
            println(io, t, ",", tm[j] / 60, ",",
                        MaxNodeDisplacement(sid; timestep = t), ",",
                        InverseCircularity(sid;  timestep = t))
        end
    end
    println("  -> ", path)
end

# ══════════════ CLAIM 2 — arbitrary / dynamic geometries ══════════════
# BLOCKED: generate_boundary_shape() is still commented out in custom.cpp and
# a/b/amp/freq are not exposed as <user_parameters>, so PCMM has no lever to
# vary the initial shape. Expose them, recompile, then flip this to true.
RUN_CLAIM2 = false

if RUN_CLAIM2
    println("Claim 2: geometry-agnostic confinement")
    GEOMETRIES = [(300.0, 300.0,  0.0, 0, "circle"),
                  (330.0, 240.0,  0.0, 0, "ellipse"),
                  (300.0, 300.0, 25.0, 4, "lobed4"),
                  (300.0, 300.0, 15.0, 8, "lobed8")]
    dv_settle = DiscreteVariation(configPath("max_time"), 1440)   # 1 day: just settle
    for (a, b, amp, freq, label) in GEOMETRIES
        out = run(inputs, dv_settle, FULL_MECHANICS..., NO_PROLIFERATION,
                  DiscreteVariation(configPath("user_parameters", "boundary_a"),    a),
                  DiscreteVariation(configPath("user_parameters", "boundary_b"),    b),
                  DiscreteVariation(configPath("user_parameters", "boundary_amp"),  amp),
                  DiscreteVariation(configPath("user_parameters", "boundary_freq"), freq);
                  n_replicates = 1, force_recompile = force_recompile)
        sid = only(simulationIDs(out))
        save_snapshot(sid, 0,      "claim2_$(label)_initial")
        save_snapshot(sid, :final, "claim2_$(label)_final")
    end
else
    println("Claim 2: skipped (RUN_CLAIM2 = false — boundary shape not yet a user_parameter)")
end

# ══════════════ CLAIM 3 — proliferation deforms the membrane inward ══════════════
# Proliferation left at the config default, so the seeded clone divides.
println("Claim 3: proliferation-driven inward deformation")

out3  = run(inputs, dv_max_time, FULL_MECHANICS...;
            n_replicates = 1, force_recompile = force_recompile)
sim3  = only(simulationIDs(out3))

MILESTONES = [125, 150, 175, 200]      # 125 = starting ring size (number_EP_cells)

for n in MILESTONES
    idx, actual = frame_at(sim3, n)
    save_snapshot(sim3, idx, "claim3_n$(actual)")
end

# Fig 3D: how the boundary changes as the clone grows.
open(joinpath(FIG_DIR, "claim3_deformation.csv"), "w") do io
    println(io, "cell_count,inverse_circularity,max_node_displacement")
    for n in MILESTONES
        _, actual = frame_at(sim3, n)
        ts = boundary_timestep_at_cell_count(sim3, n)
        println(io, actual, ",",
                    InverseCircularity(sim3;   timestep = ts), ",",
                    MaxNodeDisplacement(sim3;  timestep = ts))
    end
end
println("  -> ", joinpath(FIG_DIR, "claim3_deformation.csv"))

println("\nDone. Snapshots + CSV in $(FIG_DIR)/")
