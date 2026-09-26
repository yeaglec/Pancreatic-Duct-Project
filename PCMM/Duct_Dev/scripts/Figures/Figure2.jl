using PhysiCellModelManager, Printf

# Fresh `julia` process has no persisted project state (unlike an interactive
# session where you'd already have called this once) -- must initialize explicitly.
initializeModelManager(joinpath(@__DIR__, "..", ".."))

# TODO: rename/extend once each shape is registered as its own PCMM project
# in custom.cpp + PhysiCell_settings.xml (a/b/amp/freq exposed as user_parameters).
GEOMETRIES = [
    "DuctDev_Circle",
    "DuctDev_Star",
    # add more here as they're registered (ellipse, lobed, etc.)
]

dvs = [
    DiscreteVariation(configPath("max_time"), 7200),                         # TODO: confirm this is enough to settle per shape
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"),  1),
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"),       1),
    DiscreteVariation(configPath("user_parameters", "is_restore_lin"),      1),
    DiscreteVariation(configPath("user_parameters", "is_lumenal_pressure"), 1),
    DiscreteVariation(configPath("CAF", "cycle", "rate", 0),        0.0),     # no proliferation — isolate geometry effect alone
    DiscreteVariation(configPath("Epithelial", "cycle", "rate", 0), 0.0),
]

for proj in GEOMETRIES
    inputs = InputFolders(proj, proj; rulesets_collection = proj)   # TODO: confirm config/custom_code/rulesets share one folder name per shape

    out = run(inputs, dvs...; n_replicates = 1, force_recompile = true)
    sim_id = only(simulationIDs(out))
    println(proj, " -> sim_id = ", sim_id)

    outdir = mkpath(joinpath(@__DIR__, "figure2", proj))
    snapshot_path(idx) = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output",
                                   idx === :final ? "final.svg" : @sprintf("snapshot%08d.svg", idx))

    cp(snapshot_path(0),      joinpath(outdir, "initial.svg"); force = true)
    cp(snapshot_path(:final), joinpath(outdir, "final.svg");   force = true)
    println("  -> ", outdir)

    # TODO: containment/penetration check — no metric exists for this yet.
    # Need something like min(EP-cell distance to nearest BM node) per frame,
    # or a penetration-event counter, before this figure can report a number
    # instead of just "looks fine in the SVG."
end