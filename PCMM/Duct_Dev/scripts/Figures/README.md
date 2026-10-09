# Figure Genration Notes

## (9-27) Add Basement Membrane Plotting to PhsyiCell

PCMM's copy of PhysiCell (`PCMM/Duct_Dev/PhysiCell/modules/PhysiCell_pathology.*`) has an optional hook,
`SVG_overlay_function`, that `SVG_plot` calls to draw extra shapes on top of the cells. It is `NULL` by default.

That PhysiCell folder is its own git repo (ignored here), so the hook is saved as `physicell_svg_overlay.patch`.
On a fresh PCMM checkout, reapply it from the repo root with:
`git -C PCMM/Duct_Dev/PhysiCell apply ../scripts/Figures/physicell_svg_overlay.patch`

To draw the membrane (red line) in the SVGs of a newly imported project:

1. Make sure `draw_membrane_SVG` is in the project's `custom_modules/Utils.cpp` (declared in `Utils.h`).
2. At the end of `setup_tissue()` in `custom.cpp`, set:
   ```cpp
   SVG_overlay_function = draw_membrane_SVG;
   ```
3. Delete the old `project` executable in `data/inputs/custom_codes/<project>/` so PCMM recompiles.

Note: base PhysiCell (used by Studio) does not have this hook, so keep that line out of builds that run there.

## (9-30) Parallel division for cancer cells

`parallel_cell_division` was set on `cell_defaults` in `create_cell_types()` (`custom.cpp`). In PhysiCell the first
cell definition in the XML (Epithelial, ID 0) *is* `cell_defaults`; every other type (including CAF, our cancer cells)
is a copy made when the XML is read, before that line runs. So only Epithelial cells divided parallel to the membrane.
It is now set on each definition (Epithelial and CAF), in both `custom.cpp` copies (PCMM and `PhysiCell/user_projects`).

Effect (Figure 1, one cancer cell): both daughters stay on the membrane and push the epithelial neighbour off
(~30 min); cancer cells start stacking off the membrane after 3–4 divisions (~7–8 h) and form a second layer by ~2 days.
Any figure with cancer cells changes, so old simulations in the database are stale (the Runner resets it).

## (9-30) Perturbation parameters for Figure 2

New user parameters (both XMLs, all default to 0 / off), used in `Geometry.cpp` when the ring is placed:

- `num_removed_cells` — ring cells that are not placed (a gap)
- `num_displaced_cells` — ring cells started `displaced_distance` (15 µm) further inward, off the membrane

Both pick cells in the same order as cancer cells (around ring cell 0, at 3 o'clock).

## (9-30) Healthy-duct equilibrium for Figure 2 (in progress)

With epithelial proliferation on (Figure 2 only; the XML Epithelial rate stays 0), the healthy duct did not stay at
equilibrium: extra cells were pushed off the membrane, kept dividing and formed clumps, and a 5-cell gap never closed.
Cell pressure from the outputs: all cells start at p = 0 (a burst of divisions in the first hour), ring cells sit at
p ≈ 1.7, gap-edge cells only drop to p ≈ 0.7–0.9, against the Epithelial rule's half-max of 0.5.

Changes and tests (`Model_Evaluation/Contact_Inhibition_Tests/sweep.jl`, 3 days, Epithelial rate 0.00144):

- **Membrane contact required to divide.** `phenotype_function` (`custom.cpp`, both copies; it runs after the rules)
  sets Epithelial cycle entry to 0 when the gap between the cell surface and the membrane is ≥
  `membrane_interaction_length` (10 µm), the same range as the cell–membrane adhesion.
- **Pressure rule** (Epithelial cycle entry, half_max / hill_power), cells at 0 h → 72 h:

  | rule | contact check | control | gap |
  |---|---|---|---|
  | 0.5 / 4 (current) | no | 125 → 144 | 120 → 137, clumps |
  | 1.0 / 10 | no | 125 → 172 | 120 → 149, closes in ~24 h, then clumps |
  | 1.2 / 10 | no | 125 → 250 | 120 → 222, runaway |
  | 0.5 / 4 (current) | yes | 125 → 139 | 120 → 133, closes by ~48 h as a monolayer |
  | 1.0 / 10 | yes | 125 → 151 | 120 → 149, closes in ~24 h, then clumps |

Note: PCMM reuses a stored simulation whenever the parameters match, even if the model code changed. After a model
code change, reset the database (the Runner does) before rerunning a figure or test on its own.

Follow-up (Figure 2 runs, 7 days, Epithelial rate 0.00144, contact check on; tracing each division in the control):
every division after the first hour came from a cell *on* the membrane at p ≈ 0.7–1.9. The Hill rule never fully
switches off — at p ≈ 1.7 the base rule (0.5 / 4) leaves ~0.7% of the rate, ≈ 2 divisions/day over 125 cells; each
extra daughter is squeezed off the membrane and stops dividing there (the clumps).

| rule | control | gap | displaced (15 µm) |
|---|---|---|---|
| 0.5 / 4 | 125 → 147, clumps after ~4 d | closes by ~48 h | re-attaches by ~24 h |
| 0.5 / 8 | 125 → 129 (burst in hour 1, then flat) | does not close — edge cells (p ≈ 0.7–0.9) stay inhibited | re-attaches by ~24 h |

Displaced 20 µm behaves like 15 µm; 25 µm is beyond the adhesion reach (the cells never re-attach). Final choice below.

## (9-30) Neighbour-based contact inhibition (Duct_Dev_Figs only, off by default)

**Why.** Pressure can't tell a gap-edge cell (one neighbour, p ≈ 0.8) from a loosely packed ring cell (p ≈ 0.7–0.8),
so no pressure rule both stops growth in an intact ring and lets a gap refill (see the table above).

**What.** In `custom.cpp` (`data/inputs/custom_codes/Duct_Dev_Figs` only — not `PhysiCell/user_projects`):

- `touching_neighbors(pCell)` counts cells whose centers are closer than the sum of the two radii.
- In `phenotype_function`, when `neighbor_contact_inhibition = 1`, an Epithelial cell on the membrane gets its full
  base cycle rate if it has a free side (< 2 touching neighbours) and 0 otherwise. This overrides the pressure rule.
  Cells off the membrane still never divide (the membrane-contact check above always applies).

**Toggle.** `neighbor_contact_inhibition` (user parameter, Duct_Dev_Figs XML only): `0` = off (default, pressure rule
as before), `1` = on. Only Figure 2 turns it on for now; decide later whether it belongs in the general model.
Other figures are unaffected either way (their Epithelial cycle rate is 0).

Results with it on (Figure 2, 7 days, Epithelial rate 0.00144):

- **Control:** exactly 125 cells for 7 days, a clean monolayer, no startup burst.
- **Gap:** the edge cells divide, but daughters pile up inward instead of filling the gap (120 → 132, then stops).
  Cause: `parallel_cell_division` always places the daughter counter-clockwise along the membrane (the order of the
  membrane points), so one edge divides into the gap and the other into its own neighbour. **Noted, not changed:**
  a possible fix is to place the daughter toward the free side when the toggle is on — revisit with Dr. Bergman.
  (This also means cancer-cell daughters in Figure 1 always go counter-clockwise.)
- **Displaced:** at 15 µm the cells sit 11–16 µm from the membrane for all 7 days — beyond the
  10 µm adhesion reach — held in place by their ring neighbours, which count them as touching and so never divide.
  At 7 µm (within reach) 4 of 5 re-attach within ~6 h but the offset is barely visible.
  **Kept (proof of concept):** 25 µm off (≈27 µm surface gap) with `membrane_interaction_length` raised to 30 µm for
  that row only — the cells are pulled halfway back within 1 h, 4 of 5 are flush by ~48 h, the last by day 7.

## (10-2) Duct metrics for Figure 3 (and the Figure 4 sensitivity analysis)

**Why.** Figures 1–2 show behaviour by eye; the sensitivity analysis needs each run reduced to a few numbers. Each metric
uses only the membrane outline and the cells at one time point (so it can later run on H&E) and is dimensionless.

**What.** `scripts/Metrics.jl`, one short function per metric, no package dependencies:

- **Inverse circularity** `P² / 4πA`: overall deformation (1 = circle).
- **Max indent depth** `d / R_eff`: deepest membrane point inside the equal-area circle (0 if none).
- **Off-membrane fraction**: share of cells whose surface is ≥ 10 µm (the adhesion range) from the membrane.
- **Lumen filling**: total cell area (π r²) / area inside the membrane.

`common.jl` loads a frame and returns all four: `frame_metrics(sim_id, frame)`. The old `MembraneShape.jl` /
`ParameterOptimization.jl` are untouched (the calibration scripts still use them). Dropped from the figure set: roughness,
lobe amplitude, indent extent (hard to read, redundant), area/perimeter (repeats IC), and max node displacement
(compares node *i* across frames, but membrane nodes are added during the run).

**Figure 3 runs.** Ellipse / star / circle, each as equilibrium (cancer division off) and growth (default), 7 days.
Off-membrane fraction (27–36%) and lumen filling (0.10–0.14 → 0.14–0.24) separate growth from equilibrium in every
shape; IC and indent mostly record the starting shape, with a slow rise under growth. Without division, crowding at
tight curves pushes a few cells off on day 1 (ellipse 4, star 7, circle 0).
