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

- `Scenarios/` — one file per test scenario, listing only its changes from the baseline (like the figures).
- Runner — runs every scenario for baseline and candidate.
- Checks — regression (switch off = identical output), health numbers, metrics over time.
- Report — one HTML page per comparison.

Note: PCMM reuses stored simulations whose parameters match, even after a code change. Reset the database
after changing model code.
