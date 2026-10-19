# Figure 1: TikZ schematics (A–D) stacked above the simulated extrusion strip (E).
# Edit: ZOOM (µm window), SELECTED (frames to show, offsets from the division), and the stack_rows layout.
# Output: Figure1.svg, plus explore.svg (every saved frame, zoomed) for choosing SELECTED.
include(joinpath(@__DIR__, "..", "style.jl"))

outdir = joinpath(@__DIR__, "Outputs")
explore = joinpath(outdir, "explore")

ZOOM = zoom((260, 370), (-45, 65); bar = 20)            # µm window around the cancer cell (3 o'clock), shifted up for the second layer
EXPLORE_ZOOM = zoom((200, 400), (-100, 100); bar = 50)   # wider window for the exploration sheet

# Frames are saved as k<offset>.svg, offset in frames from the first division (k+000)
offset(f) = parse(Int, basename(f)[2:end-4])
t_div = snapshot_minutes(joinpath(explore, "k+000.svg"))
minutes_after(f) = round(Int, snapshot_minutes(f) - t_div)
files = sort(filter(endswith(".svg"), readdir(explore; join = true)); by = offset)

# EXPLORATION sheet: every saved frame, zoomed — used to pick the panels below
panels = [Panel(f; title = "$(minutes_after(f)) min", transform = EXPLORE_ZOOM) for f in files]
storyboard(panels; ncols = 8, output = joinpath(@__DIR__, "explore.svg"), overwrite = true)

# Panel E: picked from the exploration sheet — (offset in 6-min frames from the first division, stage name)
SELECTED = [(-1, "Baseline"), (0, "Division"), (15, "Extrusion"), (57, "Second layer")]
frame(k) = joinpath(explore, "k$(k < 0 ? "-" : "+")$(lpad(abs(k), 3, '0')).svg")
panels = [Panel(frame(k); title = "$name ($(replace("$(minutes_after(frame(k))) min", "-" => "−")))",
                label = i == 1 ? "(E)" : "", transform = ZOOM)
          for (i, (k, name)) in enumerate(SELECTED)]
strip = storyboard(panels; legend = LEGEND, output = nothing)

# Schematics (A–D) are hand-drawn TikZ, converted once with `pdftocairo -svg name.pdf name.svg`
schematic(name) = prefix_ids(name)(read(joinpath(@__DIR__, "Schematics", "$name.svg"), String))

fig = stack_rows([[schematic("adhesion")], [schematic("extrusion")], [schematic("membrane_mechanics")], [strip]])
write(joinpath(@__DIR__, "Figure1.svg"), fig)
println("  -> ", joinpath(@__DIR__, "Figure1.svg"))
