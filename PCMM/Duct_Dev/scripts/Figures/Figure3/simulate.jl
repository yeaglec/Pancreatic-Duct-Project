# Figure 3 (simulations): the three starting geometries, each without proliferation (equilibrium: cancer-cell division off)
# and with the full default model (growth). Metrics: Metrics.jl, via frame_metrics in common.jl.
# Edit: SHAPES (membrane_shape value => folder name), MAX_TIME. Run with PCMM_NUM_PARALLEL_SIMS=6 to run all six at once.
# Output, per run in Outputs/<shape>_<condition>/:
#   t0.svg, t<MAX_TIME>.svg   first and last snapshots
#   metrics.csv               the four metrics every hour
#   membrane.csv, cells.csv   final membrane and cells (x, y, radius, off-membrane flag), for the method schematic
include(joinpath(@__DIR__, "..", "common.jl"))

MAX_TIME = 10080   # 7 days
SHAPES = Dict(0 => "ellipse", 1 => "star", 2 => "circle")
CAF_RATE = configPath("CAF", "cycle", "rate", 0)

out = run(inputs, DiscreteVariation(configPath("max_time"), MAX_TIME),
          DiscreteVariation(configPath("user_parameters", "membrane_shape"), collect(keys(SHAPES))),
          DiscreteVariation(CAF_RATE, [0.0, 0.00144]);   # 0 = equilibrium, 0.00144 = default (growth)
          n_replicates = 1)

outdir_base = fresh_dir(joinpath(@__DIR__, "Outputs"))

for sim_id in simulationIDs(out)
    # Which run is this? Read its parameters back from the settings it ran with
    xml = read(joinpath(output_dir(sim_id), "PhysiCell_settings.xml"), String)
    shape = parse(Int, match(r"<membrane_shape[^>]*>(\d+)<", xml)[1])
    caf_rate = parse(Float64, match(r"<cell_definition name=\"CAF\".*?<rate [^>]*>([\d.eE+-]+)<"s, xml)[1])
    name = SHAPES[shape] * (caf_rate == 0 ? "_equilibrium" : "_growth")
    println(name, " -> sim_id ", sim_id)

    outdir = mkpath(joinpath(outdir_base, name))
    final = frame_at(MAX_TIME)
    save_frame(sim_id, 0, joinpath(outdir, "t0.svg"))
    save_frame(sim_id, final, joinpath(outdir, "t$(MAX_TIME).svg"))

    open(joinpath(outdir, "metrics.csv"), "w") do io
        println(io, "time,IC,indent,off_membrane,filling,cells")
        for f in 0:frame_at(60):final
            v = frame_metrics(sim_id, f)
            println(io, join((f * SVG_INTERVAL, v.IC, v.indent, v.off_membrane, v.filling, v.cells), ","))
        end
    end

    m = load_membrane(sim_id, final)
    x, y, r, _ = load_cells(sim_id, final)
    writedlm(joinpath(outdir, "membrane.csv"), m, ',')
    writedlm(joinpath(outdir, "cells.csv"), [x y r off_membrane(m, x, y, r)], ',')
end
