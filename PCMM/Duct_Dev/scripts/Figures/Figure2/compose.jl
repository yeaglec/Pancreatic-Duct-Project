# Figure 2: one row per membrane geometry, columns = time points (early -> late).
# Edit: SHAPES (rows, top to bottom; names must match simulate.jl) and the panel title.
# Output: Figure2.svg
include(joinpath(@__DIR__, "..", "style.jl"))

outdir = joinpath(@__DIR__, "Outputs")

SHAPES = ["default", "star", "circle"]

panels = Panel[]
for name in SHAPES
    # Snapshots are saved as t<minutes>.svg — order them by time
    times = sort([parse(Int, chop(f; head = 1, tail = 4)) for f in readdir(joinpath(outdir, name)) if endswith(f, ".svg")])
    for t in times
        push!(panels, Panel(joinpath(outdir, name, "t$(t).svg"); title = "$name, day $(round(t / 1440; digits = 1))"))
    end
end

montage(panels; ncols = 4, legend = LEGEND, output = joinpath(@__DIR__, "Figure2.svg"), overwrite = true)
println("  -> ", joinpath(@__DIR__, "Figure2.svg"))
