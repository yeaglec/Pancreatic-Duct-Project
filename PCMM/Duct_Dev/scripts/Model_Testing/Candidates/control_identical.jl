# Negative control: an identical copy of the baseline under a new input folder name. It must reproduce the
# baseline exactly (identical to baseline = yes for every experiment), which shows the runs are repeatable and that
# any difference a real candidate shows comes from its change. Rerun whenever the pipeline or PCMM setup changes.
(
    name = "control_identical",
    inputs = "Control_Identical_V0.0.1",   # same files as Duct_Baseline_V0.0.1
    changes = [],                          # nothing on top
    check_identical = true,                # run the regression check (optional; off when left out)
)
