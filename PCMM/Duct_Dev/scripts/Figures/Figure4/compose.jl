# Figure 4: one row per swept mechanic, columns = swept values (low -> high); the default is tagged.
# Edit: SWEEPS (rows; names must match simulate.jl), DEFAULTS (config values to tag), and the panel title.
# Output: Figure4.svg
include(joinpath(@__DIR__, "..", "style.jl"))

outdir = joinpath(@__DIR__, "Outputs")

SWEEPS = ["strain", "bending", "pressure"]
DEFAULTS = Dict("strain" => 0.1, "bending" => 4.0, "pressure" => 0.01)   # Duct_Dev_Figs config defaults

panels = Panel[]
ncols = 0
for name in SWEEPS
    # Snapshots are saved as v<value>.svg — order them by value
    values = sort([parse(Float64, chop(f; head = 1, tail = 4)) for f in readdir(joinpath(outdir, name)) if endswith(f, ".svg")])
    global ncols = max(ncols, length(values))
    for v in values
        title = "$name = $v" * (v == DEFAULTS[name] ? " (default)" : "")
        push!(panels, Panel(joinpath(outdir, name, "v$(v).svg"); title))
    end
end

montage(panels; ncols, legend = LEGEND, output = joinpath(@__DIR__, "Figure4.svg"), overwrite = true)
println("  -> ", joinpath(@__DIR__, "Figure4.svg"))
