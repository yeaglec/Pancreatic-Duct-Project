# Metric Testing

Decides which metrics to keep. `Model_Testing/` tests the **model**; this folder tests the **metrics**.
The metrics live in `../Metrics.jl` (`lesion_metrics`), shared with the figures and Model_Testing.

A metric is kept when it is **correct** (known answers on synthetic shapes) and **useful** (changes with a
parameter by more than it varies between seeds, and doesn't just repeat another metric).

Run from `PCMM/Duct_Dev/scripts/`:

```bash
cd PCMM/Duct_Dev/scripts
```

## 1. Correct: `metric_tests.jl` (seconds, no simulations)

```bash
julia --project=. Metric_Testing/metric_tests.jl
```

Synthetic cases: a circle with one and two cell layers, notches of known depth and width, a notch lined by cells,
a fold away from the lesion, a rotated shape, and the same dent at 0.5 and 5 µm node spacing.

## 2. Useful: `metric_screen.jl` (9 simulations of 7 days)

```bash
julia --project=. Metric_Testing/metric_screen.jl
```

1. One cancer cell, three settings (`baseline`, `adhesion_0.1`, `bending_2`) × seeds 0, 1, 2.
2. Every metric at the first frame where the lesion reaches 5, 10 and 20 cells (lesion size, not time: a slide has
   no clock). Samples with cells outside the membrane are dropped.
3. **Signal / noise** = difference of seed means ÷ seed standard deviation. Useful if ≥ 2.
4. **Redundancy** = rank correlation between two metrics over all samples. Redundant if |r| ≥ 0.9.

Results go to `Outputs/screen/` (`samples.csv`, `signal_noise.csv`, `redundancy.csv`) and a printed summary.
Settings are at the top of the script. Simulations are stored, so a rerun after editing a metric takes seconds.
`MAX_TIME=60` in front of the command checks that everything runs.

## Adding a metric

1. Add it to `lesion_metrics` in `../Metrics.jl`.
2. Add a case with a hand-worked answer to `metric_tests.jl`.
3. Add its name to `METRICS` in `metric_screen.jl` and rerun.
4. Record the decision below.

## Metrics

All are snapshot only (measurable on H&E), local to the lesion, in cell diameters (16.8 µm), and computed on the
membrane resampled every 2 µm. Lesion = neoplastic cells (`CAF` in the config); lesion window = membrane within one
cell diameter of a lesion cell.

| Metric | Meaning | On H&E |
|---|---|---|
| `IC` | Inverse circularity of the whole duct (the one global metric) | Duct outline |
| `fold_depth` | Deepest lesion membrane point below the convex hull | Depth of an endophytic fold |
| `fold_aspect` | fold_depth ÷ width of the fold's mouth | Flat vs papillary (old PanIN-1A vs 1B) |
| `layer_thickness` | Furthest lesion cell from the membrane, in cell layers | Stratification |
| `fold_share` | fold_depth ÷ (fold_depth + extra layers): 1 = folds, 0 = stacks | Fold vs stratify |

## Decisions (screen of 2026-10-25, sims 52–60)

| Metric | Best signal / noise | Decision |
|---|---|---|
| `IC` | 9.1 (adhesion, 10 cells) | Keep: global context |
| `fold_depth` | 16.9 | Keep |
| `fold_aspect` | 14.4 | Keep, although rank r = 0.96 with fold_depth here: these folds deepen without narrowing. Needed later to constrain fold shape |
| `layer_thickness` | 4.0 | Keep |
| `fold_share` | 12.6 | Keep: summary of fold vs layering (0.07 at baseline vs 0.84 at adhesion 0.1, at 10 cells) |
| `monolayer_fraction` | 3.1 | Dropped: redundant with layer_thickness (r = −0.92) |

Halving bending (4 → 2) is barely visible: only IC (2.4) and fold_aspect (2.2) pass 2, both at 20 cells.
Bending matters little in this range before a breach. Adhesion 0.1 breaches before 20 lesion cells, so it is compared at 5 and 10 cells only.

Retired earlier: roughness, lobe amplitude and the radial / Fourier profiles in `MembraneShape.jl` (node-spacing or
starting-shape dependent; still used by the frozen calibration scripts). `max_indent_depth`, `off_membrane_fraction`
and `lumen_filling` stay in `Metrics.jl` only so Figures 3 and 4 reproduce.
