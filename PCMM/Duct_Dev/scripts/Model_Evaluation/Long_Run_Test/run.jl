# 28-day run of the default model at half the cancer-cell division rate: the default, and Gaussian force smoothing
# (sigma 3.3) with three restoring-force strengths — four runs in parallel, then a movie of each (needs `brew install librsvg`).
# Usage (from scripts/):  PCMM_NUM_PARALLEL_SIMS=4 julia --project=. Model_Evaluation/Long_Run_Test/run.jl
# Edit: BASELINE (shared settings), RUNS (sigma, restoring strength), SAVE_INTERVAL.
# Output: Outputs/<run>.mp4 and Outputs/<run>_final.svg (copies; the full outputs stay in data/outputs/simulations/<id>)
include(joinpath(@__DIR__, "..", "..", "Figures", "common.jl"))

MAX_TIME = parse(Int, get(ENV, "MAX_TIME", "40320"))   # 28 days (set MAX_TIME=60 in the shell for a quick check)
SAVE_INTERVAL = 30   # min between SVG / full-data saves: 1344 frames → ~56 s movie at 24 fps (6 min would be ~4.7 min)

BASELINE = [
    DiscreteVariation(configPath("max_time"), MAX_TIME),
    DiscreteVariation(configPath("SVG_save_interval"), SAVE_INTERVAL),
    DiscreteVariation(configPath("full_data_interval"), SAVE_INTERVAL),
    DiscreteVariation(configPath("CAF", "cycle", "rate", 0), 0.00072),   # half the default 0.00144; Epithelial rate stays 0
]

# Membrane settings for one run: smoothing sigma (0 = off; is_gaussian_smoothing is already 1) and restoring strength (0 = off)
membrane(sigma, restore) = [DiscreteVariation(configPath("user_parameters", "membrane_force_smoothing_sigma"), sigma),
                            DiscreteVariation(configPath("user_parameters", "is_restore_lin"), restore > 0 ? 1.0 : 0.0),
                            DiscreteVariation(configPath("user_parameters", "home_lin"), restore)]

RUNS = [   # (sigma, restoring strength home_lin)
    (0.0, 0.0),        # default
    (3.3, 0.0005),
    (3.3, 0.00025),
    (3.3, 0.0001),
]

# All four as one batch, so PCMM_NUM_PARALLEL_SIMS runs them at the same time
monads = [createTrial(inputs, BASELINE..., membrane(sigma, restore)...; n_replicates = 1) for (sigma, restore) in RUNS]
out = run(Sampling(monads; n_replicates = 1))

# Movie of one simulation's snapshots. (PCMM's makeMovie renders frames with ImageMagick, whose built-in SVG renderer
# on this Mac drops every outline, including the red membrane, so frames are rendered with rsvg-convert from librsvg.)
function make_movie(sim_id, dst; fps = 24)
    dir = output_dir(sim_id)
    svgs = filter(f -> startswith(f, "snapshot") && endswith(f, ".svg"), readdir(dir))
    pngs = [joinpath(dir, replace(f, ".svg" => ".png")) for f in svgs]
    for (svg, png) in zip(svgs, pngs)
        run(`rsvg-convert -b white $(joinpath(dir, svg)) -o $png`)
    end
    run(`ffmpeg -y -loglevel error -r $fps -i $(joinpath(dir, "snapshot%08d.png")) -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"
         -vcodec libx264 -pix_fmt yuv420p -tune animation -crf 15 $dst`)
    rm.(pngs)
end

# Which run is which? Read sigma and the restoring strength back from the settings each simulation ran with
function run_name(sim_id)
    xml = read(joinpath(output_dir(sim_id), "PhysiCell_settings.xml"), String)
    value(tag) = parse(Float64, match(Regex("<$tag[^>]*>([^<]+)<"), xml)[1])
    sigma, restore = value("membrane_force_smoothing_sigma"), value("is_restore_lin") == 1 ? value("home_lin") : 0.0
    return sigma == 0 && restore == 0 ? "default" : "sigma$(sigma)_restore$(restore)"
end

outdir = mkpath(joinpath(@__DIR__, "Outputs"))
for sim_id in simulationIDs(out)
    name = run_name(sim_id)
    println(name, " -> sim_id ", sim_id)
    make_movie(sim_id, joinpath(outdir, "$name.mp4"))
    save_frame(sim_id, :final, joinpath(outdir, "$(name)_final.svg"))
end
