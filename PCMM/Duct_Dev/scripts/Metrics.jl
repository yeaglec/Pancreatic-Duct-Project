# Metrics.jl — duct metrics from one saved frame (membrane polygon + cell positions), never model internals.
# Membrane = N×2 matrix of (x, y) nodes in order (boundary_t<frame>.csv). Cells = vectors x, y, radius.
# Pure geometry, no packages. Tests: Metric_Testing/metric_tests.jl.
#
# Figure metrics (Figures 3 and 4; unchanged so the figures reproduce):
#   inverse_circularity    IC = P² / (4πA): 1 for a circle, grows as the outline deforms
#   max_indent_depth       deepest membrane point inside the equal-area circle, / its radius
#   off_membrane_fraction  share of cells whose surface is ≥ `reach` µm from the membrane (pushed off the layer)
#   lumen_filling          Σ cell areas / area inside the membrane
#
# Lesion metrics (lesion_metrics, Step 2). Snapshot only, so measurable on H&E; local to the lesion; lengths in
# cell diameters; computed on the membrane resampled every 2 µm, so the model's node spacing doesn't matter.
#   lesion window    membrane within one cell diameter of a lesion (neoplastic) cell
#   IC               inverse circularity of the whole duct (the one global metric)
#   fold_depth       deepest lesion-window point below the convex hull, i.e. below the lid over the fold
#   fold_aspect      fold_depth / width of the fold's mouth: low = broad dent, high = narrow fold
#   layer_thickness  furthest lesion cell from the membrane, in cell layers (1 = monolayer)
#   fold_share       fold_depth / (fold_depth + extra layers): 1 = only folds, 0 = only stacks, NaN = neither yet

# ── Polygon basics ──
function area(m)
    n = size(m, 1)
    return abs(sum(m[i, 1] * m[mod1(i + 1, n), 2] - m[mod1(i + 1, n), 1] * m[i, 2] for i in 1:n)) / 2
end
perimeter(m) = sum(hypot(m[mod1(i + 1, size(m, 1)), 1] - m[i, 1], m[mod1(i + 1, size(m, 1)), 2] - m[i, 2]) for i in 1:size(m, 1))

function centroid(m)
    n, cx, cy, a = size(m, 1), 0.0, 0.0, 0.0
    for i in 1:n
        j = mod1(i + 1, n)
        c = m[i, 1] * m[j, 2] - m[j, 1] * m[i, 2]
        cx += (m[i, 1] + m[j, 1]) * c
        cy += (m[i, 2] + m[j, 2]) * c
        a += c
    end
    return cx / (3a), cy / (3a)
end

# Distance from point (x, y) to the closest membrane segment
function distance_to_membrane(m, x, y)
    d = Inf
    for i in 1:size(m, 1)
        j = mod1(i + 1, size(m, 1))
        ax, ay, bx, by = m[i, 1], m[i, 2], m[j, 1], m[j, 2]
        t = clamp(((x - ax) * (bx - ax) + (y - ay) * (by - ay)) / ((bx - ax)^2 + (by - ay)^2), 0, 1)
        d = min(d, hypot(x - ax - t * (bx - ax), y - ay - t * (by - ay)))
    end
    return d
end

# Is (x, y) inside the membrane? (ray casting)
function inside(m, x, y)
    n, c = size(m, 1), false
    for i in 1:n
        j = i == 1 ? n : i - 1
        if (m[i, 2] > y) != (m[j, 2] > y) && x < (m[j, 1] - m[i, 1]) * (y - m[i, 2]) / (m[j, 2] - m[i, 2]) + m[i, 1]
            c = !c
        end
    end
    return c
end

# ── Figure metrics ──
inverse_circularity(m) = perimeter(m)^2 / (4π * area(m))

function max_indent_depth(m)
    cx, cy = centroid(m)
    R = sqrt(area(m) / π)   # radius of the circle with the same area
    return max(maximum(R - hypot(m[i, 1] - cx, m[i, 2] - cy) for i in 1:size(m, 1)) / R, 0.0)   # 0 if no point is inside (circle)
end

off_membrane(m, x, y, r; reach = 10.0) = [distance_to_membrane(m, x[k], y[k]) - r[k] >= reach for k in eachindex(x)]
off_membrane_fraction(m, x, y, r; reach = 10.0) = sum(off_membrane(m, x, y, r; reach)) / length(x)

lumen_filling(m, r) = sum(π .* r .^ 2) / area(m)

# ── Lesion metrics ──
const D_CELL = 16.8     # cell diameter (µm), the unit of length
const SPACING = 2.0     # µm between resampled membrane points
const MOUTH_EPS = 1.0   # µm: points less than this below the hull are outside the fold

# The membrane as points every `ds` µm along it
function resample(m, ds = SPACING)
    n = size(m, 1)
    s = [0.0; cumsum([hypot(m[mod1(i + 1, n), 1] - m[i, 1], m[mod1(i + 1, n), 2] - m[i, 2]) for i in 1:n])]
    k = max(round(Int, s[end] / ds), 3)
    out = zeros(k, 2)
    for q in 1:k
        t = (q - 1) * s[end] / k
        i = min(searchsortedlast(s, t), n)
        f = (t - s[i]) / (s[i + 1] - s[i])
        out[q, :] = (1 - f) * m[i, :] + f * m[mod1(i + 1, n), :]
    end
    return out
end

# Depth of every membrane point below the convex hull (0 on the hull). The hull (monotone chain) visits its vertices
# in membrane order, so each hull edge is the lid over the membrane points between its two ends.
function hull_depth(m)
    n = size(m, 1)
    turn(a, b, c) = (m[b, 1] - m[a, 1]) * (m[c, 2] - m[a, 2]) - (m[b, 2] - m[a, 2]) * (m[c, 1] - m[a, 1])
    function chain(order)
        h = Int[]
        for p in order
            while length(h) >= 2 && turn(h[end - 1], h[end], p) <= 0
                pop!(h)
            end
            push!(h, p)
        end
        return h[1:end - 1]
    end
    order = sortperm(collect(zip(m[:, 1], m[:, 2])))
    hull = sort(unique([chain(order); chain(reverse(order))]))

    depth = zeros(n)
    for (q, a) in enumerate(hull)
        b = q < length(hull) ? hull[q + 1] : hull[1] + n   # next hull vertex, wrapping around
        lid = m[mod1(b, n), :] - m[a, :]
        for i in (a + 1):(b - 1)
            v = m[mod1(i, n), :] - m[a, :]
            depth[mod1(i, n)] = abs(lid[1] * v[2] - lid[2] * v[1]) / hypot(lid...)
        end
    end
    return depth
end

# Width of the fold around point k: distance between the first points on either side that are back at the hull
function mouth_width(m, depth, k)
    n, l, r = size(m, 1), k, k
    while depth[mod1(l, n)] > MOUTH_EPS; l -= 1; end
    while depth[mod1(r, n)] > MOUTH_EPS; r += 1; end
    return hypot((m[mod1(l, n), :] - m[mod1(r, n), :])...)
end

# All lesion metrics for one frame. lesion = Bool per cell (true = neoplastic).
function lesion_metrics(m, x, y, lesion)
    m = resample(m)
    lx, ly = x[lesion], y[lesion]
    IC = inverse_circularity(m)
    fold_depth = fold_aspect = layer_thickness = fold_share = NaN
    isempty(lx) && return (; IC, fold_depth, fold_aspect, layer_thickness, fold_share)

    depth = hull_depth(m)
    window = [i for i in 1:size(m, 1) if minimum(hypot.(m[i, 1] .- lx, m[i, 2] .- ly)) <= D_CELL]
    if !isempty(window)
        k = window[argmax(depth[window])]   # deepest point of the lesion's membrane
        fold_depth = depth[k] / D_CELL
        fold_aspect = depth[k] > MOUTH_EPS ? depth[k] / mouth_width(m, depth, k) : 0.0
    end
    # a cell resting on the membrane has its centre half a diameter from it: 1 layer
    layer_thickness = maximum(distance_to_membrane(m, lx[c], ly[c]) for c in eachindex(lx)) / D_CELL + 0.5
    extra_layers = max(layer_thickness - 1, 0.0)
    fold_depth + extra_layers > 0 && (fold_share = fold_depth / (fold_depth + extra_layers))
    return (; IC, fold_depth, fold_aspect, layer_thickness, fold_share)
end
