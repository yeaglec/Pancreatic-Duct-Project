# Usage (from scripts/):  julia --project=. Figures/Runner.jl
#
# Each figure has two steps, each run as its own julia process in its own environment:
#   FigureN/simulate.jl -> runs PhysiCell, copies snapshots to FigureN/Outputs   (scripts/Project.toml, PCMM 0.3.2)
#   FigureN/compose.jl  -> builds FigureN/FigureN.svg from FigureN/Outputs       (Figures/Project.toml, Montage)

RUN_SIMS = true   # false = only recompose figures from the snapshots already in FigureN/Outputs

FIGURES = [
    "Figure1",
    "Figure2",
    "Figure3",
    "Figure4",
]

sim_env = joinpath(@__DIR__, "..")
fig_env = @__DIR__

if RUN_SIMS
    using PhysiCellModelManager
    # WARNING: wipes all prior simulations in the database
    initializeModelManager(joinpath(@__DIR__, "..", ".."))   # -> Duct_Dev
    resetDatabase(; force_reset = true, force_continue = true)
end

logdir = mkpath(joinpath(@__DIR__, "logs"))

function run_step(script, env)
    logfile = joinpath(logdir, replace(script, "/" => "_", ".jl" => ".log"))
    println("=== running $script  (log: $logfile) ===")
    try
        open(logfile, "w") do io
            run(pipeline(`$(Base.julia_cmd()) --project=$env $(joinpath(@__DIR__, script))`, stdout = io, stderr = io))
        end
        return true
    catch
        println("  FAILED: $script -- see $logfile")
        return false
    end
end

results = Dict{String, Bool}()
for fig in FIGURES
    ok = !RUN_SIMS || run_step("$fig/simulate.jl", sim_env)
    results[fig] = ok && run_step("$fig/compose.jl", fig_env)
end

println("\n=== summary ===")
for fig in FIGURES
    println(results[fig] ? "  OK   " : "  FAIL ", fig)
end
