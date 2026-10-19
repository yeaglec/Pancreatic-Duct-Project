# Healthy ring at rest: no cancer cells, and normal epithelial cells don't divide (baseline rate 0).
# Nothing should drive the membrane, so it should stay where it starts.
(
    name = "rest_ring",
    rows = [1],                 # behavior table rows this experiment tests
    max_time = 20160,           # 14 days
    changes = [                 # only what differs from the baseline config
        DiscreteVariation(configPath("user_parameters", "num_cancer_cells"), 0),
    ],
)
