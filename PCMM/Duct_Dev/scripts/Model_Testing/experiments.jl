# Loads every experiment in Experiments/. Each experiment file returns one NamedTuple:
#   name      short id, used in folder and report names
#   rows      behavior table rows it tests (see the roadmap doc)
#   max_time  run length in minutes
#   changes   DiscreteVariations: only what differs from the baseline config
# To add an experiment, drop a new file in Experiments/; nothing else changes.
load_experiments(dir = joinpath(@__DIR__, "Experiments")) =
    [include(f) for f in sort(readdir(dir; join = true)) if endswith(f, ".jl")]
