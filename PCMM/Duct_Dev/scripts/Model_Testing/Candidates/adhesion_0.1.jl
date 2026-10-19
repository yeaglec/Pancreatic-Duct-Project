# Hypothesis: stronger cell-membrane adhesion (0.01 -> 0.1) keeps an early lesion a monolayer, with the
# membrane folding to fit the new cells into the first layer.
(
    name = "adhesion_0.1",
    inputs = "Duct_Baseline_V0.0.1",   # same code and config as the baseline: this is a parameter change only
    changes = [                       # applied on top of every experiment's own changes
        DiscreteVariation(configPath("user_parameters", "membrane_adhesion_strength"), 0.1),
    ],
)
