# One cancer cell in the default duct: the smallest lesion. Does the membrane fold before a second layer
# appears, and does it fold into the lumen?
(
    name = "one_cancer_cell",
    rows = [2, 3, 4, 5],        # behavior table rows this experiment tests
    max_time = 20160,           # 14 days
    changes = [                 # only what differs from the baseline config
        DiscreteVariation(configPath("user_parameters", "num_cancer_cells"), 1),
    ],
)
