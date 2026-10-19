# Figure 2: a healthy duct at equilibrium and its recovery from each perturbation.
# Edit: ROWS (names must match simulate.jl, with the row title), HOURS (columns), ZOOM.
# Output: Figure2.svg, plus Outputs/<row>_explore.svg (whole duct every 12 h + zoom every hour for 24 h)
include(joinpath(@__DIR__, "..", "style.jl"))

outdir = joinpath(@__DIR__, "Outputs")

ROWS = [("control", "Control"), ("gap", "Gap"), ("displaced", "Displaced")]   # details (5 cells removed; 5 cells 25 µm off) go in the caption
HOURS = [0, 24, 72, 168]
ZOOM = zoom((200, 400), (-100, 100); bar = 50)   # perturbed cells sit at 3 o'clock (around ring cell 0)

frame(row, h) = joinpath(outdir, row, "h$(h).svg")
title(f) = "$(round(snapshot_minutes(f) / 60; digits = 1)) h · $(snapshot_cells(f)) cells"

# EXPLORATION sheets — used to pick the snapshot times for the figure
for (row, _) in ROWS
    whole = [Panel(frame(row, h); title = title(frame(row, h))) for h in 0:12:168 if isfile(frame(row, h))]
    zoomed = [Panel(frame(row, h); title = title(frame(row, h)), transform = ZOOM) for h in 0:24]
    storyboard([whole; zoomed]; ncols = 8, output = joinpath(outdir, "$(row)_explore.svg"), overwrite = true)
    println("  -> ", joinpath(outdir, "$(row)_explore.svg"))
end

# FIGURE: one row per condition, columns = HOURS, zoomed on the perturbed site
panels = Panel[]
for (i, (row, name)) in enumerate(ROWS), h in HOURS
    day = "day $(h ÷ 24)"
    push!(panels, Panel(frame(row, h); title = h == first(HOURS) ? "$name, $day" : day,
                        label = h == first(HOURS) ? "($('A' + i - 1))" : "", transform = ZOOM))
end
montage(panels; ncols = length(HOURS), legend = LEGEND[[1, 3]], output = joinpath(@__DIR__, "Figure2.svg"), overwrite = true)   # no cancer cells
println("  -> ", joinpath(@__DIR__, "Figure2.svg"))
