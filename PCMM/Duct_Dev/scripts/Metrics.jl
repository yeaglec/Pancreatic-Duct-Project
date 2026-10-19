# Metrics.jl — the four duct metrics used in the paper figures (introduced in Figure 3, used by the Figure 4 sensitivity analysis).
#
#   Membrane shape (from the basement-membrane polygon only):
#     inverse_circularity   IC = P² / (4πA)                     1 for a circle, grows as the outline deforms
#     max_indent_depth      max (R_eff − r) / R_eff             deepest membrane point inside the equal-area circle
#                           (r = distance of a membrane node from the centroid, R_eff = √(A/π))
#   Cells (needs the cells too):
#     off_membrane_fraction share of cells whose surface is ≥ `reach` µm from the membrane (the model's adhesion range),
#                           i.e. cells that were pushed out of the epithelial layer
#     lumen_filling         Σ cell areas (π r²) / area inside the membrane
#
# Membrane = N×2 matrix of (x, y) nodes in order (boundary_t<frame>.csv). Cells = vectors x, y, radius.
# Pure geometry, no packages, so compose scripts can use it too. Loading a simulation's frame: Figures/common.jl (frame_metrics).
# Simplified replacement for the figure use of MembraneShape.jl / ParameterOptimization.jl (those still drive the calibration scripts).

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

# ── The four metrics ──
inverse_circularity(m) = perimeter(m)^2 / (4π * area(m))

function max_indent_depth(m)
    cx, cy = centroid(m)
    R = sqrt(area(m) / π)   # radius of the circle with the same area
    return max(maximum(R - hypot(m[i, 1] - cx, m[i, 2] - cy) for i in 1:size(m, 1)) / R, 0.0)   # 0 if no point is inside (circle)
end

off_membrane(m, x, y, r; reach = 10.0) = [distance_to_membrane(m, x[k], y[k]) - r[k] >= reach for k in eachindex(x)]
off_membrane_fraction(m, x, y, r; reach = 10.0) = sum(off_membrane(m, x, y, r; reach)) / length(x)

lumen_filling(m, r) = sum(π .* r .^ 2) / area(m)
