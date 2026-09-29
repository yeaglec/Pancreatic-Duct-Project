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
