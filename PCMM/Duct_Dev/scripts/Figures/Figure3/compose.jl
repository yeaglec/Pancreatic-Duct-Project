# Figure 3: A = equilibrium (no proliferation), B = full default model; columns = the same time points.
# Edit: ROWS (top to bottom, with their sub-figure letters; names must match simulate.jl) and the panel title.
# Output: Figure3.svg
include(joinpath(@__DIR__, "..", "style.jl"))

outdir = joinpath(@__DIR__, "Outputs")

ROWS = [("equilibrium", "A"), ("growth", "B")]

panels = Panel[]
for (name, letter) in ROWS
    # Snapshots are saved as t<minutes>.svg — order them by time
    files = sort(filter(endswith(".svg"), readdir(joinpath(outdir, name); join = true)); by = snapshot_minutes)
    for (i, f) in enumerate(files)
        title = "day $(round(snapshot_minutes(f) / 1440; digits = 1)) · $(snapshot_cells(f)) cells"
        push!(panels, Panel(f; title, label = i == 1 ? letter : ""))
    end
end

montage(panels; ncols = 3, legend = LEGEND, output = joinpath(@__DIR__, "Figure3.svg"), overwrite = true)
println("  -> ", joinpath(@__DIR__, "Figure3.svg"))
