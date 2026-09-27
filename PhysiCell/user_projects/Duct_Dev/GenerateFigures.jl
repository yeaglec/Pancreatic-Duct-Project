#=
================================================================================
GenerateFigures.jl — Publication figures for Results Claims 1–3
================================================================================

Companion to ParamOpt.jl / OptSMC.jl. Where those scripts run large sweeps and
build interactive HTML comparison dashboards (GenerateReport.jl), this script
runs the smaller, targeted simulations needed for three specific manuscript
figures and composes each into a single publication-ready PDF/PNG panel.

FIGURE MAP (see results-prototyping .tex, "PART A" claims 1–3)
  Fig 1 (Claim 1) — Model schematic + baseline mechanical equilibrium.
                    "With adhesive, spring, and repulsive forces active and
                    no proliferation, the duct holds its configuration."
  Fig 2 (Claim 2) — Geometry-agnostic adhesion/repulsion across different
                    initial BM shapes. STATUS: BLOCKED, see the big comment
                    block below — the arbitrary-boundary generator exists in
                    custom.cpp but is commented out and not exposed as a
                    user_parameter, so this section is a runnable skeleton
                    for once that's wired up, not a script you can run today.
  Fig 3 (Claim 3) — Headline result: a single proliferating clone crowds and
                    extrudes toward the lumen, and the new layer pulls the BM
                    inward into a focal, growing indentation. Zoomed crop +
                    whole-duct view + quantitative depth-vs-cell-count curve.

DEPENDENCIES
  - PhysiCellModelManager (same as your other scripts)
  - CairoMakie, FileIO   (for compositing snapshot panels into each figure —
    install once with `] add CairoMakie FileIO`)
  - rsvg-convert on PATH, to rasterize PhysiCell's native SVG snapshots to PNG
    before compositing (macOS: `brew install librsvg`). If you'd rather not
    add that dependency, see the note above `svg_to_png` below for the
    ImageMagick (`convert`) alternative — same call shape, just swap the cmd.

USAGE
  Adjust the CONFIG block, then:
      julia GenerateFigures.jl
  Figures land in ./figures/ as both .pdf (for the manuscript) and .png
  (for quick viewing / slides).

A NOTE ON MAKIE API STABILITY
  `arrows!`, `text!`'s keyword names, and `Circle` have shifted a bit across
  Makie/CairoMakie versions (e.g. some releases prefer `arrows2d!`). If
  draw_force_schematic! errors on your installed version, that function is
  the one to patch — nothing else in the script depends on it.
================================================================================
=#

using PhysiCellModelManager
using CairoMakie
using FileIO
using Printf
include("GenerateReport.jl")   # also pulls in ParameterOptimization.jl:
                                #   boundary_timesteps, evaluate_simulation,
                                #   output_index_at_cell_count,
                                #   boundary_timestep_at_cell_count
                                # and GenerateReport.jl's own:
                                #   snapshot_svg_path

setNumberOfParallelSims(10)

# ══════════════════════════════════════════════════════════════════
#                    CONFIG — edit before running
# ══════════════════════════════════════════════════════════════════

const FIG_DIR = "figures"
const RASTER_DIR = joinpath(FIG_DIR, "_raster")   # intermediate PNGs, gitignore-able
mkpath(RASTER_DIR)

config_folder = custom_code_folder = rulesets_collection_folder = "DuctDev_ParamOpt"
inputs = InputFolders(config_folder, custom_code_folder;
                      rulesets_collection=rulesets_collection_folder)

dv_max_time = DiscreteVariation(configPath("max_time"), 7200)  # 5 days; bump if milestones aren't reached
force_recompile = false

# Reference "full mechanics" trial shared by Fig 1 and Fig 3 — the model as
# currently believed to be the working configuration (mirrors OptSMC.jl's
# ref_model). Re-check these against whatever your most recent calibration
# settled on before treating figures as final.
FULL_MECHANICS = [
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1),
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 1),
    DiscreteVariation(configPath("user_parameters", "is_restore_lin"), 1),
    DiscreteVariation(configPath("user_parameters", "is_lumenal_pressure"), 1),
]

# ══════════════════════════════════════════════════════════════════
#                    SHARED HELPERS
# ══════════════════════════════════════════════════════════════════

"""
Rasterize a PhysiCell SVG snapshot to PNG so it can be composited into a
Makie figure with `image!`. Requires `rsvg-convert` on PATH.

No rsvg-convert? ImageMagick works too — swap the body for:
    run(`convert -density 300 -background white $svg_path $png_path`)
"""
function svg_to_png(svg_path::String, png_path::String; width::Int=1000)
    isfile(svg_path) || error("Missing SVG: $svg_path")
    run(`rsvg-convert -w $width -o $png_path $svg_path`)
    return png_path
end

"""
    snapshot_image(sim_id, target_cell_count; label="") -> Matrix (image)

Pick the PhysiCell snapshot nearest the given cell-count milestone (or the
final frame if `target_cell_count === nothing`), rasterize it, and load it
as an image ready for `image!`. Mirrors GenerateReport.jl's own snapshot
selection logic so figures and the HTML dashboards never disagree about
which frame a given milestone means.
"""
function snapshot_image(sim_id::Int, target_cell_count::Union{Nothing,Int}; label::String="")
    if target_cell_count === nothing
        svg = snapshot_svg_path(sim_id, :final)
        tag = "final"
    else
        res = output_index_at_cell_count(sim_id, target_cell_count)
        res === nothing && error(
            "Sim $sim_id never reached $target_cell_count cells — lower the " *
            "milestone or lengthen dv_max_time.")
        idx, actual_count = res
        svg = snapshot_svg_path(sim_id, idx)
        tag = "n$(actual_count)"
        actual_count != target_cell_count &&
            @info "Sim $sim_id: target=$target_cell_count, nearest reached=$actual_count"
    end
    isfile(svg) || error("SVG not found for sim $sim_id at $tag: $svg")
    png = joinpath(RASTER_DIR, "sim$(sim_id)_$(tag)$(isempty(label) ? "" : "_"*label).png")
    svg_to_png(svg, png)
    return rotr90(load(png))   # rotr90: FileIO loads row-major/top-left origin; Makie's image! expects the flip
end

"""Bare axis with no ticks/box, for a snapshot panel. Returns the Axis."""
function image_axis(gridpos, title::String)
    ax = Axis(gridpos; title=title, titlesize=14, aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

"""Center-crop an image to a fraction of its width/height, for a 'zoomed in
on the extrusion site' panel. (x_frac, y_frac) in (0,1] control crop size;
(x_center, y_center) in [0,1] pick where on the frame to center the crop.
This is a rough starting crop — PhysiCell's SVG doesn't hand back a
domain→pixel mapping, so nudge x_center/y_center by eye until the crop
frames the extruding clone; it'll be the same seed/config every run so you
only need to tune it once and it'll stay correct as long as the deformation
occurs in roughly the same place (proliferation is seeded at a fixed cell)."""
function center_crop(img; x_frac=0.4, y_frac=0.4, x_center=0.5, y_center=0.35)
    h, w = size(img)
    cw, ch = round(Int, w * x_frac), round(Int, h * y_frac)
    cx, cy = round(Int, w * x_center), round(Int, h * y_center)
    x0 = clamp(cx - cw ÷ 2, 1, w - cw + 1)
    y0 = clamp(cy - ch ÷ 2, 1, h - ch + 1)
    return img[y0:y0+ch-1, x0:x0+cw-1]
end

"""Metric value at each cell-count milestone, for a dose-response line plot.
Reuses evaluate_simulation so the number here always matches what a report
built with GenerateReport.jl would show for the same frame."""
function metric_vs_cellcount(sim_id::Int, milestones::Vector{Int}, metric::String)
    xs, ys = Float64[], Float64[]
    for m in milestones
        res = output_index_at_cell_count(sim_id, m)
        res === nothing && continue
        _, actual = res
        ts = boundary_timestep_at_cell_count(sim_id, m)
        d = evaluate_simulation(sim_id; timestep=ts)
        haskey(d, metric) || continue
        push!(xs, actual); push!(ys, d[metric])
    end
    return xs, ys
end

# ══════════════════════════════════════════════════════════════════
#  FIGURE 1 (Claim 1) — force schematic + baseline mechanical equilibrium
# ══════════════════════════════════════════════════════════════════
println("=== Figure 1: schematic + baseline equilibrium ===")

# A single fixed configuration doesn't need createTrial — run() takes
# `inputs` plus any number of DiscreteVariations directly (confirmed by
# PCMM's own docs: `run(inputs, dv; n_replicates=3)`); createTrial is for
# building a reusable reference that a *sweep* (DistributedVariation +
# sampler) then layers variation on top of, as in ParamOpt.jl/OptSweep.jl.
out_equilibrium = run(inputs, dv_max_time, FULL_MECHANICS...,
    DiscreteVariation(configPath("user_parameters", "proliferation_exit_rate"), 0.0);  # no divisions
    n_replicates=1, force_recompile=force_recompile
)
sim_eq = only(simulationIDs(out_equilibrium))

"""Placeholder force-schematic panel: membrane arc, one basal epithelial
cell, one second-layer (lumenal) cell, and the three named forces. This is
a sketch to get the layout and labels right — swap for a hand-refined
version (or re-draw with real BM curvature/cell radii) before submission."""
function draw_force_schematic!(ax)
    θ = range(-0.55, 0.55, length=80)
    R = 300.0
    mx = R .* sin.(θ)
    my = R .* cos.(θ) .- R
    lines!(ax, mx, my; color=:black, linewidth=3)
    text!(ax, "basement membrane"; position=(mx[end]+10, my[end]), fontsize=12, align=(:left,:center))

    basal_cell  = Point2f(0, -32)   # sits on the membrane (adhesion)
    lumenal_cell = Point2f(6, -72)  # second-layer daughter, pulled toward lumen
    for (c, r, col, lbl) in ((basal_cell, 22, :steelblue, "basal EP cell"),
                              (lumenal_cell, 20, :orange, "daughter (lumenal) cell"))
        poly!(ax, Circle(c, r); color=(col, 0.35), strokecolor=col, strokewidth=2)
        text!(ax, lbl; position=(c[1]+r+8, c[2]), fontsize=11, color=col, align=(:left,:center))
    end

    bm_node = Point2f(0, my[length(my) ÷ 2])

    # Adhesive force: membrane -> basal cell
    arrows!(ax, [bm_node[1]], [bm_node[2]], [basal_cell[1]-bm_node[1]], [basal_cell[2]-bm_node[2]];
            color=:seagreen, linewidth=3, arrowsize=16)
    text!(ax, "adhesive"; position=(bm_node .+ Point2f(-70, -6)), color=:seagreen, fontsize=13)

    # Spring-like (strain) force: membrane pulled toward the deeper daughter cell
    arrows!(ax, [bm_node[1]], [bm_node[2]], [(lumenal_cell[1]-bm_node[1])*0.9], [(lumenal_cell[2]-bm_node[2])*0.9];
            color=:purple, linewidth=3, arrowsize=16, linestyle=:dash)
    text!(ax, "spring\n(deformation)"; position=(bm_node .+ Point2f(-95, -45)), color=:purple, fontsize=13)

    # Repulsive safety-net force: pushes the daughter cell back off the membrane
    arrows!(ax, [lumenal_cell[1]], [lumenal_cell[2]], [0.0], [-28.0];
            color=:firebrick, linewidth=3, arrowsize=16)
    text!(ax, "repulsive\n(safety net)"; position=(lumenal_cell .+ Point2f(15, -30)), color=:firebrick, fontsize=13)

    xlims!(ax, -140, 220)
    ylims!(ax, -140, 30)
end

fig1 = Figure(size=(1000, 500))
ax_schem = Axis(fig1[1,1], title="A. Force schematic", titlesize=14, aspect=DataAspect())
hidedecorations!(ax_schem); hidespines!(ax_schem)
draw_force_schematic!(ax_schem)

ax_eq = image_axis(fig1[1,2], "B. Equilibrium configuration (no proliferation)")
image!(ax_eq, snapshot_image(sim_eq, nothing))

save(joinpath(FIG_DIR, "fig1_equilibrium_schematic.pdf"), fig1)
save(joinpath(FIG_DIR, "fig1_equilibrium_schematic.png"), fig1, px_per_unit=3)
println("  -> figures/fig1_equilibrium_schematic.[pdf|png]")

# ══════════════════════════════════════════════════════════════════
#  FIGURE 2 (Claim 2) — geometry-agnostic adhesion/repulsion
# ══════════════════════════════════════════════════════════════════
#
#  ** BLOCKED — do not run this section until the code change below lands **
#
#  What I found in custom_modules/custom.cpp: the current setup_tissue()
#  only spawns one test CAF cell (Test_KernelCell()); the generator for
#  non-circular boundaries is already written but commented out:
#
#      // boundary_membrane_pts = generate_boundary_shape(a, b, amp, freq);
#      // generate_boundary_cells(a, b, amp, freq, "Epithelial", ep_dis, num_ep);
#
#  and a/b/amp/freq are local constants in that commented block, not
#  user_parameters — so PCMM has no lever to vary shape between sims.
#  You confirmed the ring setup (generate_circle_cells) works in the
#  DuctDev_ParamOpt copy of the code, but this arbitrary-boundary path was
#  not confirmed working there either — treat it as needing the same fix.
#
#  Suggested fix (small, does not touch mechanics code):
#    1. In config/PhysiCell_settings.xml <user_parameters>, add:
#         <boundary_a    type="double" .../>  (semi-major axis, default = membrane_circle_radius)
#         <boundary_b    type="double" .../>  (semi-minor axis, default = membrane_circle_radius, i.e. circle)
#         <boundary_amp  type="double" .../>  (ripple amplitude, default 0)
#         <boundary_freq type="int"    .../>  (ripple frequency, default 0)
#    2. In custom.cpp's setup_tissue(), replace the commented block with:
#         double a = parameters.doubles("boundary_a");
#         double b = parameters.doubles("boundary_b");
#         double amp = parameters.doubles("boundary_amp");
#         int freq = parameters.ints("boundary_freq");
#         int num_ep = parameters.ints("number_EP_cells");
#         double ep_dis = parameters.doubles("ep_displacement");
#         boundary_membrane_pts = generate_boundary_shape(a, b, amp, freq);
#         generate_boundary_cells(a, b, amp, freq, "Epithelial", ep_dis, num_ep);
#       (keep this ABOVE the existing `initialize_level_set_duct(...)` call,
#       and remove/guard whatever currently produces the ring so the two
#       don't both run)
#    3. Recompile DuctDev_ParamOpt.
#
#  Once that's live, this is the sweep: circle (amp=0) vs. a couple of
#  amplitude/frequency combinations, each run to a fixed *time* (not cell
#  count — Claim 2 is about the healthy/pre-proliferation shape, so there's
#  no milestone to match to) and snapshotted once settled.
# ══════════════════════════════════════════════════════════════════

RUN_FIGURE_2 = false   # flip to true once the code change above is in place

if RUN_FIGURE_2
    println("=== Figure 2: geometry-agnostic adhesion/repulsion ===")

    # (a, b, amp, freq, label) — first entry is the circular control.
    GEOMETRIES = [
        (300.0, 300.0, 0.0,  0, "circle"),
        (330.0, 240.0, 0.0,  0, "ellipse"),
        (300.0, 300.0, 25.0, 4, "4-lobed"),
        (300.0, 300.0, 15.0, 8, "8-lobed (rough)"),
    ]
    dv_settle_time = DiscreteVariation(configPath("max_time"), 1440)  # 1 day: just needs to settle, no growth

    fig2 = Figure(size=(1400, 400))
    for (i, (a, b, amp, freq, label)) in enumerate(GEOMETRIES)
        out = run(inputs, dv_settle_time, FULL_MECHANICS...,
            DiscreteVariation(configPath("user_parameters", "proliferation_exit_rate"), 0.0),
            DiscreteVariation(configPath("user_parameters", "boundary_a"), a),
            DiscreteVariation(configPath("user_parameters", "boundary_b"), b),
            DiscreteVariation(configPath("user_parameters", "boundary_amp"), amp),
            DiscreteVariation(configPath("user_parameters", "boundary_freq"), freq);
            n_replicates=1, force_recompile=force_recompile
        )
        sid = only(simulationIDs(out))
        ax = image_axis(fig2[1, i], "$(('A':'Z')[i]). $label")
        image!(ax, snapshot_image(sid, nothing; label=label))
    end

    save(joinpath(FIG_DIR, "fig2_geometry_agnostic.pdf"), fig2)
    save(joinpath(FIG_DIR, "fig2_geometry_agnostic.png"), fig2, px_per_unit=3)
    println("  -> figures/fig2_geometry_agnostic.[pdf|png]")
else
    println("=== Figure 2 skipped (RUN_FIGURE_2 = false — see code-change note above) ===")
end

# ══════════════════════════════════════════════════════════════════
#  FIGURE 3 (Claim 3) — aberrant proliferation drives membrane deformation
# ══════════════════════════════════════════════════════════════════
println("=== Figure 3: proliferation-driven focal deformation (headline) ===")

# Matched cell-count milestones: an early one (crowding/extrusion just
# starting) and a later one (clear focal indentation). number_EP_cells=125
# is the starting ring size in the current config — adjust if yours differs.
const MILESTONE_EARLY = 140
const MILESTONE_LATE  = 200
const ALL_MILESTONES  = [125, 135, 140, 150, 165, 180, 200]  # for the dose-response curve

# This run deliberately leaves proliferation_exit_rate at its config default
# (nonzero) — that default is what drives the single seeded cell's clone.
# If proliferation onset in your model is instead a discrete on/off switch
# elsewhere (rather than the pressure-modulated exit rate in
# cell_rules.csv), add that DiscreteVariation here instead.
out_proliferating = run(inputs, dv_max_time, FULL_MECHANICS...;
    n_replicates=1, force_recompile=force_recompile
)
sim_prolif = only(simulationIDs(out_proliferating))

img_early = snapshot_image(sim_prolif, MILESTONE_EARLY; label="early")
img_late  = snapshot_image(sim_prolif, MILESTONE_LATE;  label="late")
img_zoom  = center_crop(img_early)   # tune center_crop's x_center/y_center to frame the extrusion site

fig3 = Figure(size=(1300, 650))
ax_zoom  = image_axis(fig3[1,1], "A. Extrusion & crowding (n=$MILESTONE_EARLY, zoom)")
image!(ax_zoom, img_zoom)
ax_early = image_axis(fig3[1,2], "B. Whole duct (n=$MILESTONE_EARLY)")
image!(ax_early, img_early)
ax_late  = image_axis(fig3[1,3], "C. Focal deformation (n=$MILESTONE_LATE)")
image!(ax_late, img_late)

xs, ys = metric_vs_cellcount(sim_prolif, ALL_MILESTONES, "max_indent_depth")
ax_metric = Axis(fig3[2, 1:3], xlabel="cell count", ylabel="max indentation depth (µm)",
                  title="D. Deformation grows with continued proliferation")
scatterlines!(ax_metric, xs, ys; color=:firebrick, markersize=12, linewidth=2)

save(joinpath(FIG_DIR, "fig3_focal_deformation.pdf"), fig3)
save(joinpath(FIG_DIR, "fig3_focal_deformation.png"), fig3, px_per_unit=3)
println("  -> figures/fig3_focal_deformation.[pdf|png]")

println("\nDone. Figures 1 and 3 generated in '$(FIG_DIR)/'; Figure 2 is a skeleton pending the boundary-generator wiring described above.")
