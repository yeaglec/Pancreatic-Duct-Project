# Model Testing Pipeline

Tests whether a change to the model (a new mechanic, a parameter, a code cleanup) improves the basement
membrane without breaking what already works. Every test is a **comparison**: the current **baseline** model vs a
**candidate**, run on the same set of **experiments**, scored by the same **checks**, written up as one HTML report.

Plan and behavior table: [Membrane Accuracy Roadmap](https://claude.ai/code/artifact/547a1a8a-ea04-4e8e-baf3-6bf344db6621).

```
Experiments/*.jl ─┐
                  ├─► run_comparison.jl ─► simulations ─► checks.jl ─► report.jl ─► Outputs/<candidate>/report.html
Candidates/X.jl ──┘      (baseline + candidate, 1 sim each per experiment)
```

---

## Quick start

All commands run from `PCMM/Duct_Dev/scripts/`.

```bash
cd PCMM/Duct_Dev/scripts
```

Run a comparison (here the `adhesion_0.1` candidate), then open the report:

```bash
julia --project=. Model_Testing/run_comparison.jl adhesion_0.1
```

```bash
open Model_Testing/Outputs/adhesion_0.1/report.html
```

Quick check that everything runs (every experiment cut to 60 min):

```bash
MAX_TIME=60 julia --project=. Model_Testing/run_comparison.jl adhesion_0.1
```

Rebuild only the report from saved results (no simulations):

```bash
julia --project=. Model_Testing/report.jl adhesion_0.1
```

Baseline simulations are stored and reused, so each new comparison only simulates the candidate.

---

## What's in this folder

| Path | What it is | Edit it when |
|---|---|---|
| `Experiments/` | One file per experiment: a test situation (e.g. one cancer cell). Lists only its changes from the baseline | Adding a test situation |
| `Candidates/` | One file per candidate: which input folder to run, plus parameter changes | Testing anything new |
| `run_comparison.jl` | Runs everything, then checks and report. Settings at the top | Changing the baseline version or run settings |
| `checks.jl` | Computes every number from the saved outputs. Thresholds at the top | Adding or tuning a check |
| `report.jl` | Builds the HTML report | Changing the report |
| `experiments.jl` | Loads `Experiments/` | Rarely |
| `run_baseline.jl` | One 7-day reference run of the baseline | Rarely |
| `Outputs/<candidate>/` | `report.html`, `summary.csv`, `checks/*.csv`, `simulations.csv` (git-ignored) | Never: generated |

Model inputs live in `PCMM/Duct_Dev/data/inputs/` as three folders with the same name:
`configs/<name>/`, `custom_codes/<name>/`, `rulesets_collections/<name>/`.

---

## Recipes

### 1. Test a parameter change

No code change, so the candidate reuses the baseline's input folder. Create `Candidates/bending_2.jl`:

```julia
# Hypothesis: softer bending (4 -> 2) lets the membrane fold more tightly around the lesion.
(
    name = "bending_2",
    inputs = "Duct_Baseline_V0.0.1",     # same code and config as the baseline
    changes = [                         # applied on top of every experiment
        DiscreteVariation(configPath("user_parameters", "membrane_bending_constant"), 2.0),
    ],
)
```

```bash
julia --project=. Model_Testing/run_comparison.jl bending_2
```

The comment at the top becomes the hypothesis in the report. `changes` can list several parameters.

### 2. Test a new mechanic

A mechanic changes code, so it gets **its own input folder** and sits **behind a switch** (0 = off).

**Step 1: copy the baseline inputs to a new name.** Copies the source files only, not build files.

```bash
cd ../data/inputs
OLD=Duct_Baseline_V0.0.1
NEW=Candidate_WallContact
for k in configs custom_codes rulesets_collections; do
  rsync -a --exclude='project' --exclude='compilation.*' --exclude='*.db' --exclude='*_variations' \
        --exclude='physicell_commit_hash.txt' $k/$OLD/ $k/$NEW/
done
cd ../../scripts
```

**Step 2: add the switch and its parameters** to `data/inputs/configs/Candidate_WallContact/PhysiCell_settings.xml`,
inside `<user_parameters>`:

```xml
<is_wall_contact type="double" units="none">0</is_wall_contact>
<wall_contact_distance type="double" units="micron">8.4</wall_contact_distance>
```

**Step 3: write the mechanic** in `data/inputs/custom_codes/Candidate_WallContact/custom_modules/Membrane.cpp`.
Call it inside `update_basement_membrane_deformation`, before the node positions are updated:

```cpp
if (parameters.doubles("is_wall_contact") == 1) membrane_wall_contact(node_forces);
```

```cpp
// Keep non-neighbouring membrane walls at least wall_contact_distance apart
void membrane_wall_contact(std::vector<std::pair<double,double>>& node_forces)
{
    double d_min = parameters.doubles("wall_contact_distance");
    // ... add forces to node_forces ...
}
```

and declare it in `Membrane.h`:

```cpp
void membrane_wall_contact(std::vector<std::pair<double,double>>& node_forces);
```

**Step 4: prove the switch-off version changes nothing.** Create `Candidates/wall_contact_off.jl`:

```julia
# Wall contact code added but switched off: must be identical to the baseline.
(
    name = "wall_contact_off",
    inputs = "Candidate_WallContact",
    changes = [],
    check_identical = true,              # adds the "identical to baseline?" table to the report
)
```

```bash
julia --project=. Model_Testing/run_comparison.jl wall_contact_off
```

The report must say **identical: yes** for every experiment. If not, the new code changed something it shouldn't have.

**Step 5: test the mechanic.** Create `Candidates/wall_contact.jl`:

```julia
# Hypothesis: wall contact stops folds from pinching shut, so the membrane stays continuous (rows 4-5).
(
    name = "wall_contact",
    inputs = "Candidate_WallContact",
    changes = [
        DiscreteVariation(configPath("user_parameters", "is_wall_contact"), 1.0),
        DiscreteVariation(configPath("user_parameters", "wall_contact_distance"), 8.4),
    ],
)
```

```bash
julia --project=. Model_Testing/run_comparison.jl wall_contact
```

To try a few values of the new parameter, make one candidate file per value (e.g. `wall_contact_d17.jl`).

> If you edit the code in `Candidate_WallContact` again, either use a new folder name or reset the database
> (below). PCMM reuses a stored simulation when the folder name and parameters match, even if the code changed.

### 3. Test a code cleanup (no intended change)

Same as recipe 2, steps 1 and 4: copy the inputs, make the cleanup, and run a candidate with
`check_identical = true` and no `changes`. **Identical: yes** means the cleanup changed nothing.

### 4. Add an experiment

Create a file in `Experiments/`. It runs automatically in every comparison from then on. Example
`Experiments/three_cancer_cells.jl`:

```julia
# Three neighbouring cancer cells: a wider lesion.
(
    name = "three_cancer_cells",
    rows = [2, 3, 4, 5],                 # behavior table rows it tests
    max_time = 20160,                    # 14 days, in minutes
    changes = [                          # only what differs from the baseline config
        DiscreteVariation(configPath("user_parameters", "num_cancer_cells"), 3),
    ],
)
```

The baseline runs once for the new experiment (stored afterwards); each comparison then adds its candidate run.

### 5. Approve a candidate: new baseline version

When a candidate is accepted, it becomes the next baseline. Old baselines are kept and never edited.

1. Copy the candidate's inputs to the new version name with the step-1 command from recipe 2
   (`OLD=Candidate_WallContact`, `NEW=Duct_Baseline_V0.0.2`).
2. In the new config, set the switch's default to on (e.g. `<is_wall_contact ...>1</is_wall_contact>`).
3. Set `BASELINE = "Duct_Baseline_V0.0.2"` at the top of `run_comparison.jl` and `run_baseline.jl`.
4. Add a row to the version table below.
5. Copy the code and config into `PhysiCell/user_projects/Duct_Dev`, keeping its base-PhysiCell differences (below).

### 6. Reset the database

Only needed after editing code inside an existing input folder. It deletes all stored simulations, so the
next comparison reruns the baseline too. Answer `y` at the prompt.

```bash
julia --project=. -e 'using PhysiCellModelManager; initializeModelManager(".."); resetDatabase()'
```

---

## Experiments

| Experiment | What happens | Rows tested |
|---|---|---|
| `rest_ring` | Healthy ring, no cancer cells, no division. The membrane should stay still | 1 |
| `one_cancer_cell` | One dividing cancer cell in the default duct | 2, 3, 4, 5 |

---

## Reading the report

1. **Header**: the hypothesis and the settings that differ between baseline and candidate.
2. **Summary table**: one row pair per experiment. Green / red = candidate better / worse.
3. **Identical to baseline?**: only when the candidate sets `check_identical = true`.
4. **Snapshots**: baseline above candidate at matching moments. A red border means the frame is after a breach.
5. **Lesion metrics at matched lesion size**: each metric when the lesion first reaches 5, 10 and 20 cells, before
   any breach. This is the quantitative comparison; it uses lesion size, not time, as on a slide.
6. **Charts**: baseline (grey) and candidate (blue) on the same axes. Dashed lines mark each run's first breach;
   values after it include escaped cells and are not meaningful.

## What the checks measure

All computed from saved outputs only (membrane file + cell positions). Per frame unless noted.

| Check | Meaning | Concern when |
|---|---|---|
| Cells outside | Cells whose centre is outside the membrane | Any for 2 consecutive frames (a **breach**) |
| Closest wall gap | Closest approach of two membrane points ≥ 40 µm apart along the membrane | Below a cell radius (8.4 µm): walls touching |
| Displacement | Distance of the membrane from its starting shape | Large at rest (`rest_ring`) |
| Stretch | Perimeter, and largest edge, relative to the start | Large: the BM should bend, not stretch |
| Curvature radius | Smallest radius of a circle through 3 membrane points 8.4 µm apart | Far below a cell size: spikes |
| Order of events *(per sim)* | Does the membrane fold (moves ≥ 8.4 µm since the first division) before a second-layer cell appears? | Fold first = no |
| Lesion metrics | IC, fold depth, fold aspect, layer thickness, fold share ([definitions](../Metric_Testing/README.md)) | Context, not pass/fail |

Thresholds are constants at the top of `checks.jl` (`R_CELL`, `ARC_EXCLUDE`, `CURV_SPAN`, `REACH`, `BREACH_FRAMES`,
`LESION_SIZES`).

## Settings

| Setting | File | Default |
|---|---|---|
| `BASELINE` | `run_comparison.jl`, `run_baseline.jl` | `Duct_Baseline_V0.0.1` |
| `NUM_PARALLEL_SIMS` | `run_comparison.jl` | 4 |
| `SAVE_INTERVAL` (min between saved frames) | `run_comparison.jl` | 30 |
| `BREACH_FRAMES` (frames outside before a breach counts) | `checks.jl` | 2 |

---

## Baseline versions

| Version | Date | Source | Change |
|---|---|---|---|
| V0.0.1 | 2026-10-19 | `Duct_Dev_Figs` (Dev-Cam `ccba464`) | Starting point. Unchanged copy, minus two uncompiled scratch files (`Models copy.cpp`, `custom copy.cpp`). |

`Control_Identical_V0.0.1` is an identical copy of V0.0.1, used by `Candidates/control_identical.jl` (the
negative control). It passed on 2026-10-19: two separate runs of the same code are identical. Rerun it only if the
pipeline or the PCMM setup changes.

`PhysiCell/user_projects/Duct_Dev` is kept in sync with the current baseline for Studio runs. It differs only where
base PhysiCell's API differs from PCMM's fork: `main.cpp`, two renamed calls in `custom.cpp`
(`setup_cell_rules`, `load_cells_from_pugixml`), and the SVG membrane overlay line, which is commented out.
