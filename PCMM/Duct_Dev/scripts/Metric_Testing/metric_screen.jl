# Is each lesion metric useful? Runs a few settings × seeds and measures every metric at matched lesion sizes.
#   signal / noise  difference between two settings' seed means ÷ the seed standard deviation (useful if ≥ SNR_MIN)
#   redundancy      rank correlation between two metrics across all samples (redundant if |r| ≥ R_MAX)
# Usage (from scripts/):  julia --project=. Metric_Testing/metric_screen.jl     (MAX_TIME=60 for a quick check)
# Output: Outputs/screen/{samples, signal_noise, redundancy}.csv and a printed summary.
# PCMM stores the simulations, so after editing a metric a rerun only recomputes the numbers.
using PhysiCellModelManager, Statistics, DelimitedFiles

initializeModelManager(joinpath(@__DIR__, "..", ".."))   # -> Duct_Dev
include(joinpath(@__DIR__, "..", "Metrics.jl"))

# ── Settings ──
INPUTS = "Duct_Baseline_V0.0.1"
MAX_TIME = parse(Int, get(ENV, "MAX_TIME", "10080"))   # 7 days
SAVE_INTERVAL = 30
SEEDS = [0, 1, 2]
LESION_SIZES = [5, 10, 20]   # lesion cell counts to compare at (a slide has no clock)
LESION_TYPE = "CAF"          # cell definition of the neoplastic cells
METRICS = (:IC, :fold_depth, :fold_aspect, :layer_thickness, :fold_share)
SNR_MIN, R_MAX = 2.0, 0.9
setNumberOfParallelSims(4)

user(p, v) = DiscreteVariation(configPath("user_parameters", p), v)
SETTINGS = ["baseline" => [],
            "adhesion_0.1" => [user("membrane_adhesion_strength", 0.1)],
            "bending_2" => [user("membrane_bending_constant", 2.0)]]
COMPARISONS = [("baseline", "adhesion_0.1"), ("baseline", "bending_2")]
COMMON = [user("num_cancer_cells", 1),   # one dividing cancer cell
          DiscreteVariation(configPath("max_time"), MAX_TIME),
          DiscreteVariation(configPath("SVG_save_interval"), SAVE_INTERVAL),
          DiscreteVariation(configPath("full_data_interval"), SAVE_INTERVAL)]

# ── Run: one simulation per setting and seed ──
inputs = InputFolders(INPUTS, INPUTS; rulesets_collection = INPUTS)
runs = [(setting = name, seed = seed,
         monad = createTrial(inputs, COMMON..., changes..., DiscreteVariation(["options", "random_seed"], seed); n_replicates = 1))
        for (name, changes) in SETTINGS for seed in SEEDS]
run(Trial([Sampling([r.monad for r in runs])]))

# ── Measure: every metric at the first frame where the lesion reaches each size ──
output_dir(id) = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$id", "output")
function load_frame(id, f)
    m = readdlm(joinpath(output_dir(id), "boundary_t$(f).csv"), ',', Float64)[:, 1:2]
    c = PhysiCellSnapshot(id, f; include_cells = true).cells
    c = c[(.!c.dead) .& (c.cell_type_name .!= "BM"), :]
    return m, c.position_1, c.position_2, c.cell_type_name .== LESION_TYPE
end

samples = []
for r in runs
    id = only(simulationIDs(r.monad))
    todo = copy(LESION_SIZES)
    f = 0
    while !isempty(todo) && isfile(joinpath(output_dir(id), "boundary_t$(f).csv"))
        m, x, y, lesion = load_frame(id, f)
        while !isempty(todo) && count(lesion) >= todo[1]
            outside = count(k -> !inside(m, x[k], y[k]), eachindex(x))   # escaped cells: sample not used
            push!(samples, merge((setting = r.setting, seed = r.seed, sim_id = id, lesion_size = popfirst!(todo),
                                  day = f * SAVE_INTERVAL / 1440, outside = outside), lesion_metrics(m, x, y, lesion)))
        end
        f += 1
    end
end
usable = filter(s -> s.outside == 0, samples)
isempty(usable) && (println("no usable samples: runs too short for the lesion to grow (fine for a MAX_TIME check)"); exit())

# ── Score ──
values_of(k; setting = nothing, size = nothing) =
    [s[k] for s in usable if isfinite(s[k]) && setting in (nothing, s.setting) && size in (nothing, s.lesion_size)]
ranks(v) = [count(<(a), v) + (count(==(a), v) + 1) / 2 for a in v]   # ties share their mean rank

signal = []
for (a, b) in COMPARISONS, n in LESION_SIZES, k in METRICS
    va, vb = values_of(k; setting = a, size = n), values_of(k; setting = b, size = n)
    (length(va) < 2 || length(vb) < 2) && continue   # need two seeds of each to estimate the noise
    noise = sqrt((var(va) + var(vb)) / 2)
    push!(signal, (metric = k, comparison = "$a vs $b", lesion_size = n, mean_a = mean(va), mean_b = mean(vb),
                   noise = noise, snr = abs(mean(vb) - mean(va)) / noise))
end

redundancy = []
for (i, p) in enumerate(METRICS), q in METRICS[(i + 1):end]
    pairs = [(s[p], s[q]) for s in usable if isfinite(s[p]) && isfinite(s[q])]
    r = cor(ranks(first.(pairs)), ranks(last.(pairs)))
    push!(redundancy, (metric_a = p, metric_b = q, n = length(pairs), rank_r = r, redundant = abs(r) >= R_MAX))
end

# ── Write and print ──
function write_csv(path, rows)
    open(path, "w") do io
        println(io, join(keys(first(rows)), ","))
        foreach(r -> println(io, join((v isa AbstractFloat ? round(v; sigdigits = 4) : v for v in r), ",")), rows)
    end
end
outdir = mkpath(joinpath(@__DIR__, "Outputs", "screen"))
write_csv(joinpath(outdir, "samples.csv"), samples)
write_csv(joinpath(outdir, "signal_noise.csv"), signal)
write_csv(joinpath(outdir, "redundancy.csv"), redundancy)

println("\n$(length(usable)) of $(length(samples)) samples usable (the rest had cells outside the membrane)")
println("\nBest signal / noise (useful if ≥ $SNR_MIN):")
for k in METRICS
    best = argmax(s -> s.snr, filter(s -> s.metric == k, signal))
    println("  ", rpad(k, 17), lpad(round(best.snr; digits = 1), 6), "   ", best.comparison, " at ", best.lesion_size,
            " cells   ", best.snr >= SNR_MIN ? "useful" : "weak")
end
println("\nRedundant pairs (|rank r| ≥ $R_MAX):")
foreach(r -> println("  ", r.metric_a, " ~ ", r.metric_b, "   ", round(r.rank_r; digits = 2)), filter(r -> r.redundant, redundancy))
println("\nfiles: ", outdir)
