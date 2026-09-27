using PhysiCellModelManager

# WARNING: wipes prior output — see note above on resetDatabase vs deleteAllSimulations
resetDatabase(force_reset = true)

SCRIPTS = [
    "Figure1.jl",
    "Figure3.jl",
    "Figure4.jl",
    # "Figure2.jl",  # not ready — needs custom.cpp geometry work first
]

logdir = mkpath(joinpath(@__DIR__, "logs"))
results = Dict{String, Bool}()

for script in SCRIPTS
    logfile = joinpath(logdir, replace(script, ".jl" => ".log"))
    println("=== running $script  (log: $logfile) ===")
    ok = try
        open(logfile, "w") do io
            run(pipeline(`julia $(joinpath(@__DIR__, script))`, stdout = io, stderr = io))
        end
        true
    catch
        println("  FAILED: $script -- see $logfile")
        false
    end
    results[script] = ok
end

println("\n=== summary ===")
for (script, ok) in results
    println(ok ? "  OK   " : "  FAIL ", script)
end