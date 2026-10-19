# Figure 3: how we measure the duct.
#   A  the four metrics, drawn on one simulated duct (EXAMPLE's last frame)
#   B  first and last snapshot of each starting geometry (growth), with their metric values
#   C  the metrics every hour for 7 days: growth (solid) vs equilibrium (dashed), one colour per geometry
# Edit: SHAPES (folder prefix, name, colour), EXAMPLE, METRICS (order, titles, y labels).
# Output: Figure3.svg
include(joinpath(@__DIR__, "..", "style.jl"))
include(joinpath(@__DIR__, "..", "..", "Metrics.jl"))
using Printf

outdir = joinpath(@__DIR__, "Outputs")

SHAPES = [("ellipse", "Ellipse", "#0072B2"), ("star", "Star", "#E69F00"), ("circle", "Circle", "#009E73")]
EXAMPLE = "ellipse_growth"
METRICS = [(:IC,           "Inverse circularity", "IC"),
           (:indent,       "Max indent depth",    "d / R_eff"),
           (:off_membrane, "Off-membrane cells",  "fraction of cells"),
           (:filling,      "Lumen filling",       "fraction of lumen area")]
DAYS = 7

# Rows A and C are drawn at the same scale as B's storyboard (6 × 300 px panels), so all text ends up the same size
W, TITLE, BODY, SMALL = 450, 22, 18, 15
REACH = 10   # µm, the model's cell–membrane adhesion range (membrane_interaction_length)

# ── Reading the simulate.jl outputs ──
parse_rows(lines) = reduce(vcat, [parse.(Float64, split(l, ","))' for l in lines])
readmatrix(f) = parse_rows(readlines(f))
function readmetrics(name)                                        # metrics.csv -> Dict(column => values)
    lines = readlines(joinpath(outdir, name, "metrics.csv"))
    vals = parse_rows(lines[2:end])
    d = Dict(Symbol(c) => vals[:, i] for (i, c) in enumerate(split(lines[1], ",")))
    d[:indent] = max.(d[:indent], 0)   # runs saved before Metrics.jl clamped this (a circle gives −0.0003)
    return d
end

txt(x, y, s; size = BODY, anchor = "middle", weight = "normal", fill = "black", extra = "") =
    """<text x="$x" y="$y" text-anchor="$anchor" font-family="Arial" font-size="$size" font-weight="$weight" fill="$fill" $extra>$s</text>"""
heading(title, label) = txt(W / 2, TITLE, title; size = TITLE, weight = "bold") *
                        (isempty(label) ? "" : txt(0, TITLE, label; size = TITLE, anchor = "start", weight = "bold"))

# ── A: the metrics, drawn on one duct ──
m = readmatrix(joinpath(outdir, EXAMPLE, "membrane.csv"))
cells = readmatrix(joinpath(outdir, EXAMPLE, "cells.csv"))      # x, y, radius, off-membrane (1/0)
cx, cy = centroid(m)
R = sqrt(area(m) / π)
deep = argmin([hypot(m[i, 1] - cx, m[i, 2] - cy) for i in 1:size(m, 1)])   # deepest membrane node
DRAW = 360                                                                   # drawing box (px)
span = 2 * maximum(max(abs(m[i, 1] - cx), abs(m[i, 2] - cy)) for i in 1:size(m, 1)) + 30
sc = DRAW / span
px(x, y) = (W / 2 + sc * (x - cx), TITLE + 16 + DRAW / 2 - sc * (y - cy))  # µm -> px, centred on the centroid
pts(m) = join(("$(round(a; digits = 1)),$(round(b; digits = 1))" for (a, b) in (px(m[i, 1], m[i, 2]) for i in 1:size(m, 1))), " ")
membrane(; fill = "none") = """<polygon points="$(pts(m))" fill="$fill" stroke="red" stroke-width="2.5"/>"""
ref_circle() = (c = px(cx, cy); """<circle cx="$(c[1])" cy="$(c[2])" r="$(sc * R)" fill="none" stroke="black" stroke-width="1.5" stroke-dasharray="7 5"/>""")
cell_dots(color) = join(("""<circle cx="$(px(c[1], c[2])[1])" cy="$(px(c[1], c[2])[2])" r="$(sc * c[3])" fill="$(color(c[4] == 1))" stroke="#555" stroke-width="0.5"/>"""
                         for c in eachrow(cells)), "\n")
function schematic(body, title, formula, note; label = "")
    y = TITLE + 16 + DRAW
    return """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$(y + 70)">
    $(heading(title, label))
    $body
    $(txt(W / 2, y + 32, formula))
    $(txt(W / 2, y + 58, note; size = SMALL, fill = "#555"))
    </svg>"""
end

c0 = px(cx, cy)
pdeep = px(m[deep, 1], m[deep, 2])
u = (m[deep, 1] - cx, m[deep, 2] - cy) ./ hypot(m[deep, 1] - cx, m[deep, 2] - cy)
pring = px(cx + R * u[1], cy + R * u[2])                         # where the ray to the deepest node meets the circle
prad = px(cx + R * u[2], cy - R * u[1])                          # a radius at right angles, labelled R_eff
dot(p, r = 4) = """<circle cx="$(p[1])" cy="$(p[2])" r="$r" fill="black"/>"""
line(a, b; w = 2, dash = "") = """<line x1="$(a[1])" y1="$(a[2])" x2="$(b[1])" y2="$(b[2])" stroke="black" stroke-width="$w" $(isempty(dash) ? "" : "stroke-dasharray=\"$dash\"")/>"""
mid(a, b) = ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2)

n_off = Int(sum(cells[:, 4]))
row_A = [
    schematic(membrane(fill = "#fde7e7") * ref_circle() * txt(c0[1], c0[2] + 7, "A"; size = TITLE, weight = "bold"),
              "Inverse circularity", "IC = P² / 4πA", "P perimeter, A area · 1 = circle (dashed)"; label = "(A)"),
    schematic(membrane() * ref_circle() * dot(c0) * line(c0, prad; w = 1.5, dash = "4 3") *
              txt(mid(c0, prad)..., "R_eff"; size = SMALL, extra = "dy=\"-6\"") *
              line(pring, pdeep; w = 4) * dot(pdeep, 5) * txt(mid(pring, pdeep)..., "d"; size = TITLE, weight = "bold", extra = "dx=\"14\""),
              "Max indent depth", "max indent = d / R_eff", "deepest point inside the equal-area circle"),
    schematic(membrane() * cell_dots(off -> off ? "#1a1a1a" : "#e6e6e6"),
              "Off-membrane cells", "N_off / N  ( = $n_off / $(size(cells, 1)) here)", "black: cell surface ≥ $REACH µm from the membrane"),
    schematic(membrane(fill = "#d6e8f7") * cell_dots(_ -> "#8c8c8c"),
              "Lumen filling", "Σ cell area / A", "cells (grey) over the area inside the membrane (blue)"),
]

# ── B: first and last snapshot per geometry (growth), metric values underneath ──
CROP = 690                                                    # µm window; stops above PhysiCell's own scale bar (y < −350)
ZOOM = zoom((-345, 345), (-340, 350); bar = 200, font = 38)
BAND = 44   # metric text size, in snapshot px (the crop is scaled down ~4×)
function with_metrics(v)
    return function (svg)
        inner = ZOOM(svg)
        return """<svg xmlns="http://www.w3.org/2000/svg" width="$CROP" height="$(CROP + 2.8BAND)">
        $inner
        $(txt(CROP / 2, CROP + 1.15BAND, @sprintf("IC %.2f  ·  indent %.2f", v.IC, v.indent); size = BAND))
        $(txt(CROP / 2, CROP + 2.4BAND, @sprintf("off-membrane %.0f%%  ·  filling %.2f", 100v.off_membrane, v.filling); size = BAND))
        </svg>"""
    end
end

panels = Panel[]
for (j, (shape, name, _)) in enumerate(SHAPES)
    mt = readmetrics("$(shape)_growth")
    for (k, (t, day)) in enumerate(((0, 0), (DAYS * 1440, DAYS)))
        i = findfirst(==(t), mt[:time])
        v = (IC = mt[:IC][i], indent = mt[:indent][i], off_membrane = mt[:off_membrane][i], filling = mt[:filling][i])
        push!(panels, Panel(joinpath(outdir, "$(shape)_growth", "t$t.svg");
                            title = k == 1 ? "$name, day 0" : "day $day · $(Int(mt[:cells][i])) cells",
                            label = j == 1 && k == 1 ? "(B)" : "", transform = with_metrics(v)))
    end
end
row_B = storyboard(panels; legend = LEGEND, output = nothing)

# ── C: metrics over time ──
function nice_ticks(lo, hi; n = 4)
    step = 10.0^floor(log10((hi - lo) / n))
    for k in (1, 2, 5, 10)
        (hi - lo) / (k * step) <= n && (step *= k; break)
    end
    return collect((floor(lo / step) * step):step:(ceil(hi / step) * step + step / 2))
end

H = 340
L, Rm, T, B = 78, 14, TITLE + 22, 58                              # plot margins
series = Dict((s, c) => readmetrics("$(s)_$(c)") for (s, _, _) in SHAPES, c in ("growth", "equilibrium"))
function plot(key, title, ylabel; label = "")
    ys = vcat((v[key] for v in values(series))...)
    yt = nice_ticks(minimum(ys), maximum(ys))
    X(d) = L + (W - L - Rm) * d / DAYS
    Y(v) = H - B - (H - B - T) * (v - yt[1]) / (yt[end] - yt[1])
    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H">""", heading(title, label))
    for v in yt
        println(io, """<line x1="$L" y1="$(Y(v))" x2="$(W - Rm)" y2="$(Y(v))" stroke="#e5e5e5"/>""",
                txt(L - 8, Y(v) + 5, rstrip(rstrip(@sprintf("%.3f", v), '0'), '.'); size = SMALL, anchor = "end"))
    end
    for d in 0:DAYS
        println(io, txt(X(d), H - B + 22, d; size = SMALL))
    end
    println(io, """<rect x="$L" y="$T" width="$(W - L - Rm)" height="$(H - B - T)" fill="none" stroke="black"/>""",
            txt(L + (W - L - Rm) / 2, H - 8, "time (days)"; size = SMALL),
            txt(0, 0, ylabel; size = SMALL, extra = "transform=\"translate(18,$((T + H - B) / 2)) rotate(-90)\""))
    for (s, _, color) in SHAPES, c in ("equilibrium", "growth")
        v = series[(s, c)]
        p = join(("$(round(X(t / 1440); digits = 1)),$(round(Y(y); digits = 1))" for (t, y) in zip(v[:time], v[key])), " ")
        println(io, """<polyline points="$p" fill="none" stroke="$color" stroke-width="2.5" $(c == "equilibrium" ? "stroke-dasharray=\"7 5\"" : "")/>""")
    end
    println(io, "</svg>")
    return String(take!(io))
end
row_C = [plot(k, t, y; label = i == 1 ? "(C)" : "") for (i, (k, t, y)) in enumerate(METRICS)]

# Key for C, as wide as row C so it is scaled the same
KW = 4W + 3 * 12
key = IOBuffer()
x = KW / 2 - 470
for (_, name, color) in SHAPES
    print(key, """<line x1="$x" y1="20" x2="$(x + 40)" y2="20" stroke="$color" stroke-width="4"/>""", txt(x + 50, 26, name; anchor = "start"))
    global x += 150
end
x += 30
print(key, """<line x1="$x" y1="20" x2="$(x + 40)" y2="20" stroke="black" stroke-width="2.5"/>""", txt(x + 50, 26, "growth (cancer cells divide)"; anchor = "start"))
x += 300
print(key, """<line x1="$x" y1="20" x2="$(x + 40)" y2="20" stroke="black" stroke-width="2.5" stroke-dasharray="7 5"/>""", txt(x + 50, 26, "equilibrium (no division)"; anchor = "start"))
key_svg = """<svg xmlns="http://www.w3.org/2000/svg" width="$KW" height="40">$(String(take!(key)))</svg>"""

fig = stack_rows([row_A, [row_B], row_C, [key_svg]])
write(joinpath(@__DIR__, "Figure3.svg"), fig)
println("  -> ", joinpath(@__DIR__, "Figure3.svg"))
