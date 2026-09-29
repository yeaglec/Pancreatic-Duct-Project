# Figure 1: TikZ schematics (A–D) stacked above the simulated extrusion strip (E).
# Edit: ZOOM (µm window), SELECTED (frames to show, offsets from the division), and the stack_rows layout.
# Output: Figure1.svg, plus explore.svg (every saved frame, zoomed) for choosing SELECTED.
include(joinpath(@__DIR__, "..", "style.jl"))

outdir = joinpath(@__DIR__, "Outputs")
explore = joinpath(outdir, "explore")

ZOOM = zoom((260, 370), (-55, 55); bar = 20)   # µm window around the cancer cell (3 o'clock)

# Frames are saved as k<offset>.svg, offset in frames from the first division (k+00)
offset(f) = parse(Int, basename(f)[2:4])
t_div = snapshot_minutes(joinpath(explore, "k+00.svg"))
minutes_after(f) = round(Int, snapshot_minutes(f) - t_div)
files = sort(filter(endswith(".svg"), readdir(explore; join = true)); by = offset)

# EXPLORATION sheet: every saved frame, zoomed — used to pick the extrusion frames below
panels = [Panel(f; title = "$(minutes_after(f)) min", transform = ZOOM) for f in files]
storyboard(panels; ncols = 8, output = joinpath(@__DIR__, "explore.svg"), overwrite = true)

# Panel E: extrusion after the first division, picked from the exploration sheet (offsets in 6-min frames)
SELECTED = [-1, 0, 4, 10]   # before, division, extruding, extruded
picked = [f for f in files if offset(f) in SELECTED]
panels = [Panel(f; title = replace("$(minutes_after(f)) min", "-" => "−"), label = i == 1 ? "(E)" : "", transform = ZOOM)
          for (i, f) in enumerate(picked)]
strip = storyboard(panels; legend = LEGEND, output = nothing)

# Schematics (A–D) are hand-drawn TikZ, converted once with `pdftocairo -svg name.pdf name.svg`
schematic(name) = prefix_ids(name)(read(joinpath(@__DIR__, "Schematics", "$name.svg"), String))

fig = stack_rows([[schematic("adhesion")], [schematic("extrusion")], [schematic("membrane_mechanics")], [strip]])
write(joinpath(@__DIR__, "Figure1.svg"), fig)
println("  -> ", joinpath(@__DIR__, "Figure1.svg"))
