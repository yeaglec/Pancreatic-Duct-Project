# Model Testing Pipeline

Tests whether a model change (a new mechanic) improves the basement membrane without breaking what already works.
Plan and behavior table: [Membrane Accuracy Roadmap](https://claude.ai/code/artifact/547a1a8a-ea04-4e8e-baf3-6bf344db6621).

## Baseline and versioning

Every comparison is **baseline vs candidate**. Both are PCMM inputs in `data/inputs/` (config, custom code, rules).

- **Baseline:** `Duct_Baseline_V<major>.<minor>.<patch>`. The current one is `Duct_Baseline_V0.0.1`.
- **Candidate:** a new input folder holding the baseline plus one change, with the new mechanic behind a
  parameter switch (0 = off, so the candidate with the switch off reproduces the baseline exactly).
- When a candidate is approved, it becomes the next baseline version. Old baselines are kept, never edited.

| Version | Date | Source | Change |
|---|---|---|---|
| V0.0.1 | 2026-10-19 | `Duct_Dev_Figs` (Dev-Cam `ccba464`) | Starting point. Unchanged copy, minus two uncompiled scratch files (`Models copy.cpp`, `custom copy.cpp`). |

`PhysiCell/user_projects/Duct_Dev` is kept in sync with the current baseline for Studio runs. It differs only where
base PhysiCell's API differs from PCMM's fork: `main.cpp`, two renamed calls in `custom.cpp`
(`setup_cell_rules`, `load_cells_from_pugixml`), and the SVG membrane overlay line, which is commented out.

## Layout (built step by step)

- `Experiments/` — one file per experiment, listing only its changes from the baseline (like the figures).
- `Candidates/` — one file per candidate: its input folder and any parameter changes on top of each experiment.
- `run_comparison.jl <candidate>` — runs every experiment for baseline and candidate (1 simulation each, 30-min
  saves) and writes `Outputs/<candidate>/simulations.csv`. `MAX_TIME=60` in the shell makes a quick check.
- `checks.jl` — computed from saved outputs only, called at the end of `run_comparison.jl`:
  per-frame health numbers (cells outside, wall gap, longest edge, displacement), numerical checks (stretch,
  curvature), Metrics.jl numbers → `Outputs/<candidate>/checks/<experiment>_<version>.csv`; per-simulation
  order of events and first breach (2 consecutive frames) → `Outputs/<candidate>/summary.csv`;
  `same_output(sim_a, sim_b)` for the regression check.
- Report — one HTML page per comparison. *(next)*

A candidate that changes code needs its own input folder name: PCMM reuses stored simulations whose inputs and
parameters match, even if the code inside a folder changed.

Note: PCMM reuses stored simulations whose parameters match, even after a code change. Reset the database
after changing model code.
