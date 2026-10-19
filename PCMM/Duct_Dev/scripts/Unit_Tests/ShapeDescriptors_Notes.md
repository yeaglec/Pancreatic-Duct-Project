# Basement-Membrane Shape Descriptors — Notes

*Companion to `MembraneShape.jl`. Explains what summary statistics we added to the calibration
pipeline, why we chose them, how they are defined, and why they behave the way they do.*

---

## 1. The problem these statistics solve

The calibration pipeline compares a **simulated** duct to a **target** duct by reducing each membrane
shape to a handful of numbers (summary statistics), then minimizing the distance between the two sets
of numbers. The quality of the whole calibration is therefore capped by the quality of these numbers.

The previous statistics were **`IC`** (inverse circularity) and **`area_frac_change`**. Two problems:

1. **They are global.** They describe the whole ring at once. Our key experiment — one cell
   proliferating and pushing a *local* dent into the membrane — is a *focal* deformation. A single
   global number can't tell "one deep local dent" apart from "mild deformation spread all the way
   around," yet those are biologically very different (a nascent lesion vs. uniform dilation).
2. **They under-constrain the fit.** In an unweighted mean-squared-error distance, `IC` (deviations
   ~0.1) dominated `area_frac_change` (~0.05), so calibration was effectively driven by *one* number
   against *three* parameters. That is non-identifiable: many parameter sets reproduce one number, so
   the "best fit" is a smear across the prior rather than a point.

**Design constraint — measurability.** We are building for eventual **H&E histology data**, where we
have a *single* traced lumen contour, not a simulation with tracked nodes or a known initial state.
So every statistic must be computable from **one closed contour**, using only quantities intrinsic to
that contour. This rules out per-node displacement, "membrane strain," or anything needing the
initial shape. It is the reason we chose the representation below.

---

## 2. The unifying representation: the radial profile r(θ)

Instead of treating the membrane as a list of 1000 indexed nodes, we describe it as a **radius as a
function of angle**, measured from the shape's own center. Concretely, for a contour with vertices
$(x_i, y_i)$:

1. **Centroid** $(c_x, c_y)$ — the area-weighted polygon centroid (not the vertex average, which
   would be biased by uneven node spacing). Code: `polygon_centroid`.
2. For each vertex, its **angle** $\theta_i = \operatorname{atan2}(y_i - c_y,\; x_i - c_x)$ and
   **radius** $r_i = \sqrt{(x_i-c_x)^2 + (y_i-c_y)^2}$.
3. **Resample** onto a uniform grid of $M = 256$ angles by sorting on $\theta$ and interpolating
   (periodically, so the seam at $\pm\pi$ joins). This gives a clean function $r(\theta)$.
   Code: `radial_profile`.

Why this representation is the right one:

- **Robust to remeshing.** The model inserts/removes membrane nodes, so node *index* i is not stable
  across frames. But the *contour* is. Describing the shape by $r(\theta)$ (a geometric object)
  sidesteps node identity entirely.
- **H&E-measurable.** Tracing a lumen boundary from an image and sampling its radius vs. angle is
  exactly what a pathologist/segmentation tool can produce. Same code runs on sim and image.
- **Naturally separates local from global.** A dent is a localized dip in $r(\theta)$; overall
  dilation is a change in its mean; multi-lobing is periodic structure. Different biological behaviors
  live in different features of the same curve.

**Scale normalizer.** We define the **effective radius** $R_{\text{eff}} = \sqrt{A/\pi}$ — the radius
of a circle with the same area $A$. Every length is divided by $R_{\text{eff}}$, so all descriptors
are **dimensionless** and comparable between a 270 µm duct and a 50 µm duct (and between sim and
image). This is computable from the single contour (it only needs the enclosed area).

---

## 3. The descriptors

All are returned by `shape_descriptors(pts)`. Let $r(\theta)$ be the resampled profile and
$\tilde r(\theta) = r(\theta)/R_{\text{eff}}$ the normalized profile.

### 3.1 IC — inverse circularity *(global, kept)*

$$\mathrm{IC} = \frac{P^2}{4\pi A}$$

where $P$ = perimeter, $A$ = area. **Why it works:** the isoperimetric inequality says $P^2 \ge 4\pi A$
for *any* closed curve, with equality **only** for a circle. So $\mathrm{IC} = 1$ for a circle and
$>1$ for anything else — a clean, unit-free "how non-circular is this." It is translation-, rotation-,
and scale-invariant. **Limitation:** it says *how much* deviation there is, not *where* or *what kind*.
That's what the next three add.

### 3.2 max_indent_depth — deepest local inward dip *(local)*

Define the pointwise inward depth relative to the equal-area circle:
$$d(\theta) = \frac{R_{\text{eff}} - r(\theta)}{R_{\text{eff}}} = 1 - \tilde r(\theta),
\qquad \text{max\_indent\_depth} = \max_\theta d(\theta).$$

**What it captures:** the single deepest inward push, as a fraction of duct size. This is your
"greatest deformation," made scale-free and image-measurable. **Why it works:** the equal-area circle
is the natural neutral reference derivable from the contour alone (no initial frame needed); a focal
dent is precisely where $r$ falls well below that reference. **Behavior:** ≈ 0 for a circle; for a
Gaussian dent of depth $d_0$ it recovers ≈ $d_0/R_{\text{eff}}$ (slightly less if the dent is narrow
and under-sampled). This is the descriptor that most directly tracks the proliferation experiment.

### 3.3 indent_extent — angular width of deformation *(local)*

$$\text{indent\_extent} = \frac{\#\{\theta : d(\theta) > \delta\}}{M}, \qquad \delta = 0.05.$$

I.e. the fraction of the ring that is indented more than a 5 % threshold. **What it captures:** *how
wide* the deformation is — one narrow crypt-like dent vs. deformation spread around the ring. Paired
with `max_indent_depth` (how deep) it distinguishes **focal** from **distributed** deformation, which
`IC` alone cannot. **H&E-measurable:** it's just the arc length of the "pinched in" region.
**Caveat (see §5):** with our freq-4 starting geometry, the baseline already deviates from a circle,
so this metric has a non-zero baseline and is less sensitive to the *incremental* dent than
`max_indent_depth`.

### 3.4 lobe_amp and roughness — the shape's "frequency content"

We take the Fourier spectrum of the normalized profile $\tilde r(\theta)$ (after subtracting its
mean). Intuitively, we ask *"what periodic patterns make up this contour?"* A wavenumber $k$ means "a
pattern that repeats $k$ times around the ring": $k=2$ is an ellipse, $k=4$ is a four-lobed clover,
high $k$ is fine jaggedness. The amplitude at $k$ says how strong that pattern is. Code:
`fourier_amplitudes` (a direct DFT — we avoid an FFTW dependency; the signal is small, $M=256$).

- **lobe_amp** = the largest amplitude among **low** wavenumbers $k = 2\ldots6$. **Captures**
  organized, large-scale lobing / multi-site deformation. (Our initial geometry is a $k=4$ shape, so
  its `lobe_amp` sits around the initial amplitude 0.1.)
- **roughness** = the **RMS amplitude of high** wavenumbers $k > 10$:
  $\sqrt{\sum_{k>10} a_k^2}$. **Captures** fine, disorganized jaggedness — the fingerprint of a
  membrane going numerically unstable (nodes zig-zagging). It is a cheap early-warning / validity
  signal.

  *Why it is defined as an absolute amplitude and not a fraction of total energy:* a fraction blows up
  for a smooth contour (a perfect circle's tiny numerical noise is broadband, so a large *fraction* of
  a near-zero total sits at high $k$ — giving a meaningless "roughness = 0.3"). The absolute RMS is ~0
  for a circle and grows only with real high-frequency structure. This was a bug we caught and fixed.

---

## 4. Where we measure them: the cell-count milestone

Previously statistics were read at the **final** simulation frame. But the final frame might catch a
duct mid-collapse — not a fair, reproducible comparison point. We now evaluate at a **fixed
cell-count milestone**: the first frame at which the duct reaches *N* cells (`target_cell_count`).
This compares every parameter set at a **comparable biological state** ("given N cells crowding the
duct, where has the membrane settled?") rather than at an arbitrary wall-clock time. Code:
`boundary_timestep_at_cell_count` (shared with the report so both pick the same frame).

Practical note: the milestone must be a count the sims actually reach (ours reach ~138–193 cells, so
we use 130). If it isn't reached, the code warns and falls back to the final frame.

---

## 5. Honest caveats (say these to your advisor before they ask)

- **`indent_extent` has a non-zero baseline** because our initial membrane is a 4-lobed shape, not a
  circle. So it partly measures the *starting* geometry, not just the new dent. Measuring deviation
  from the *initial* contour would isolate the lesion better — but the initial contour doesn't exist
  in a single H&E image, so we deliberately kept the single-contour definition. The focal signal shows
  up cleanly in `max_indent_depth` regardless.
- **`area_frac_change` is the one statistic that needs an initial frame** ($ (A_t - A_0)/A_0 $). It is
  fine for simulations but will **not** transfer to a single histology image. When we move to real
  data we'll swap in a single-contour occlusion measure (e.g. lumen area ÷ convex-hull area).
- These are **shape** descriptors of the membrane only. They don't yet use cell positions (e.g. how
  many epithelial layers have formed). That's a deliberate later addition.
- Adding more/better statistics improves *identifiability*, but the calibration still needs
  **normalization + a sensibly chosen tolerance** to converge well — that's the next stage of work,
  not something the statistics alone fix.

---

## 6. Evidence they behave correctly

**Analytic ground truth** (`TestDescriptors.jl`, 34 checks): a perfect circle gives IC ≈ 1 and all
deformation descriptors ≈ 0; a $\cos(k\theta)$ ripple of amplitude $\varepsilon$ produces a spectral
peak exactly at wavenumber $k$ with amplitude ≈ $\varepsilon$; a Gaussian dent recovers its depth and
width; everything is invariant to translating and scaling the shape.

**On real simulations** (all start from the identical initial shape, so t0 values match exactly):

| sim / frame | IC | max_indent_depth | indent_extent | lobe_amp | roughness |
|---|---|---|---|---|---|
| initial (all) | 1.08 | 0.141 | 0.336 | 0.111 | 0.003 |
| sim 7 final | 1.333 | 0.247 | 0.344 | 0.091 | 0.037 |
| sim 21 final | 1.286 | **0.295** | 0.336 | 0.093 | 0.036 |
| sim 4 final (stable) | 1.084 | 0.173 | 0.324 | 0.111 | 0.005 |

The key line is **sim 21 vs sim 7**: sim 21 has a *deeper* focal dent (0.295 vs 0.247) even though its
*global* circularity is *lower* (IC 1.286 vs 1.333). A global-only calibration would rank these
backwards for "worst local deformation." That single row is the argument for why the local descriptors
were worth adding.
