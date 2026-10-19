# One self-contained HTML page per comparison: Outputs/<candidate>/report.html
# Built from the files run_comparison.jl writes (simulations.csv, summary.csv, checks/*.csv, candidate.jl) plus the
# simulations' snapshots and settings. No JavaScript; charts are inline SVG, snapshots are embedded.
# Standalone (from scripts/):  julia --project=. Model_Testing/report.jl <candidate>
using Printf, Base64, Dates

# ── Reading the CSVs ──
function read_table(path)   # Vector of Dicts, values parsed to Float64 where possible ("missing" -> missing)
    lines = readlines(path)
    header = split(lines[1], ",")
    parsev(v) = v == "missing" ? missing : v in ("true", "false") ? v == "true" : something(tryparse(Float64, v), String(v))
    return [Dict(String(h) => parsev(v) for (h, v) in zip(header, split(l, ","))) for l in lines[2:end] if !isempty(l)]
end
column(rows, name) = [r[name] for r in rows]

# ── Formatting ──
fmt(v) = ismissing(v) ? "—" : v isa Bool ? (v ? "yes" : "no") : v isa Real ? (isfinite(v) ? @sprintf("%.3g", v) : "—") : string(v)
html(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

# Which way is better for the coloured columns (others are shown without a verdict)
const BETTER = Dict("fold_first" => :true, "first_breach" => :later, "max_stretch" => :lower,
                    "max_edge_stretch" => :lower, "min_curv_radius" => :higher)
function verdict(name, b, c)   # "better" / "worse" / "" for the candidate cell
    rule = get(BETTER, name, nothing)
    (rule === nothing || isequal(b, c)) && return ""
    rule == :true && return ismissing(c) || ismissing(b) ? "" : c ? "better" : "worse"
    if rule == :later   # missing = never breached = best
        ismissing(c) && return "better"
        ismissing(b) && return "worse"
        return c > b ? "better" : "worse"
    end
    (ismissing(b) || ismissing(c) || !isfinite(b) || !isfinite(c)) && return ""
    return (rule == :lower) == (c < b) ? "better" : "worse"
end

# ── Settings diff (same structure in both files: compare line by line) ──
function settings_diff(sim_a, sim_b)
    clean(l) = strip(replace(l, r" description=\"[^\"]*\"" => ""))
    la = clean.(readlines(joinpath(output_dir(sim_a), "PhysiCell_settings.xml")))
    lb = clean.(readlines(joinpath(output_dir(sim_b), "PhysiCell_settings.xml")))
    length(la) == length(lb) || return [("(settings files differ in structure)", "", "")]
    rows = Tuple{String, String, String}[]
    for (a, b) in zip(la, lb)
        a == b && continue
        ma, mb = match(r"<([\w]+)[^>]*>([^<]*)<", a), match(r"<([\w]+)[^>]*>([^<]*)<", b)
        push!(rows, ma === nothing || mb === nothing ? (a, "", b) : (ma[1], ma[2], mb[2]))
    end
    return rows
end

# ── Snapshots ──
snapshot_file(sim_id, frame) = joinpath(output_dir(sim_id), @sprintf("snapshot%08d.svg", frame))
function snapshot_img(sim_id, frame, label, breached)
    f = snapshot_file(sim_id, frame)
    isfile(f) || return "<figure class=\"snap\"><div class=\"missing\">no snapshot</div><figcaption>$(html(label))</figcaption></figure>"
    data = base64encode(read(f))
    cls = breached ? "snap breached" : "snap"
    return "<figure class=\"$cls\"><img src=\"data:image/svg+xml;base64,$data\" alt=\"$(html(label))\"><figcaption>$(html(label))</figcaption></figure>"
end
function frame_of_day(rows, day)   # closest saved frame to a day (day may be missing)
    ismissing(day) && return nothing
    k = argmin(abs.(column(rows, "day") .- day))
    return Int(rows[k]["frame"])
end
function snapshot_row(sim_id, rows, summary)
    last_day = rows[end]["day"]
    picks = [("start", 0.0), ("first division", summary["first_division"]), ("first 2nd layer", summary["first_second_layer"]),
             ("day $(fmt(last_day / 4))", last_day / 4), ("day $(fmt(last_day / 2))", last_day / 2),
             ("day $(fmt(3last_day / 4))", 3last_day / 4), ("final", last_day)]
    breach = summary["first_breach"]
    cells = String[]
    for (label, day) in picks
        f = frame_of_day(rows, day)
        f === nothing && (push!(cells, "<figure class=\"snap\"><div class=\"missing\">never</div><figcaption>$label</figcaption></figure>"); continue)
        d = rows[f + 1]["day"]
        caption = startswith(label, "day") ? label : "$label · day $(fmt(d))"
        push!(cells, snapshot_img(sim_id, f, caption, !ismissing(breach) && d >= breach))
    end
    return join(cells)
end

# ── Charts: baseline vs candidate on shared axes ──
function chart(title, unit, b_rows, c_rows, col, b_breach, c_breach; from_zero = true)
    W, H, L, R, T, B = 340, 200, 48, 12, 28, 30
    xs = vcat(column(b_rows, "day"), column(c_rows, "day"))
    ys = filter(isfinite, Float64.(vcat(column(b_rows, col), column(c_rows, col))))
    isempty(ys) && return ""
    x0, x1 = 0.0, maximum(xs)
    y0, y1 = from_zero ? min(0.0, minimum(ys)) : minimum(ys), maximum(ys)
    y1 == y0 && (y1 = y0 + 1)
    X(x) = L + (x - x0) / (x1 - x0) * (W - L - R)
    Y(y) = H - B - (y - y0) / (y1 - y0) * (H - T - B)
    function path(rows)
        d, pen = IOBuffer(), false
        for r in rows
            v = r[col]
            if v isa Real && isfinite(v)
                print(d, pen ? "L" : "M", @sprintf("%.1f %.1f", X(r["day"]), Y(v))); pen = true
            else
                pen = false
            end
        end
        return String(take!(d))
    end
    breach_line(day, cls) = ismissing(day) ? "" :
        @sprintf("<line class=\"%s\" x1=\"%.1f\" x2=\"%.1f\" y1=\"%d\" y2=\"%d\"/>", cls, X(day), X(day), T, H - B)
    ticks = join(@sprintf("<text x=\"%d\" y=\"%.1f\" class=\"tick\" text-anchor=\"end\">%s</text><line class=\"grid\" x1=\"%d\" x2=\"%d\" y1=\"%.1f\" y2=\"%.1f\"/>",
                          L - 6, Y(v) + 4, fmt(v), L, W - R, Y(v), Y(v)) for v in range(y0, y1; length = 4))
    xt = join(@sprintf("<text x=\"%.1f\" y=\"%d\" class=\"tick\" text-anchor=\"middle\">%d</text>", X(d), H - B + 16, d)
              for d in 0:7:floor(Int, x1))
    return """<svg viewBox="0 0 $W $H" class="chart" role="img" aria-label="$(html(title))">
      <text x="$L" y="16" class="ctitle">$(html(title))$(isempty(unit) ? "" : " ($(html(unit)))")</text>$ticks$xt
      <text x="$(W - R)" y="$(H - 4)" class="tick" text-anchor="end">days</text>
      $(breach_line(b_breach, "breach-b"))$(breach_line(c_breach, "breach-c"))
      <path d="$(path(b_rows))" class="line-b"/><path d="$(path(c_rows))" class="line-c"/></svg>"""
end

const CHARTS = [("Closest wall gap", "µm", "wall_gap", true), ("Max displacement from start", "µm", "disp_max", true),
                ("Perimeter stretch", "× start", "stretch", false), ("Smallest curvature radius", "µm", "min_curv_radius", true),
                ("Off-membrane fraction", "", "off_membrane", true), ("Cells", "", "cells", true), ("Cells outside", "", "outside", true)]

const NOTES = [
    ("Cells outside", "Cells whose centre is outside the membrane. A breach is 2 consecutive frames with any; frames after it have a red border."),
    ("Closest wall gap", "Closest approach of two membrane points ≥ 40 µm apart along the membrane: two walls of a fold touching. Below one cell radius (8.4 µm) is a warning."),
    ("Displacement", "Distance of the membrane from its starting contour. At rest (row 1) it should stay near zero."),
    ("Perimeter / edge stretch", "Perimeter, and the largest edge, relative to the start. The BM should fold by bending, so large stretch is suspect. Edge stretch is only valid while remeshing is off."),
    ("Curvature radius", "Smallest radius of a circle through three membrane points 8.4 µm apart. Radii far below a cell size are spikes."),
    ("Order of events", "First fold: membrane moves ≥ 8.4 µm from where it was at the first division. First 2nd layer: off-membrane count rises above its value at the first division. Row 2 passes if the fold comes first."),
    ("Off-membrane fraction", "Share of cells whose surface is ≥ 10 µm (the adhesion range) from the membrane."),
]

const CSS = """
:root { --ink:#1d1d1f; --quiet:#6e6e73; --rule:#e5e5ea; --bg:#ffffff; --panel:#f6f6f8; --base:#8e8e93; --cand:#2463eb;
        --good:#1a7f37; --bad:#c62828; }
@media (prefers-color-scheme: dark) { :root { --ink:#f2f2f7; --quiet:#a1a1a6; --rule:#3a3a3c; --bg:#1c1c1e; --panel:#2c2c2e;
        --base:#a1a1a6; --cand:#6ea0ff; --good:#4cc26e; --bad:#ff6b6b; } }
body { font: 14px/1.45 -apple-system, system-ui, sans-serif; color: var(--ink); background: var(--bg); margin: 0; }
main { max-width: 1180px; margin: 0 auto; padding: 24px 16px 64px; }
h1 { font-size: 22px; margin: 0 0 4px; } h2 { font-size: 17px; margin: 32px 0 8px; border-bottom: 1px solid var(--rule); padding-bottom: 4px; }
h3 { font-size: 15px; margin: 20px 0 6px; } .quiet { color: var(--quiet); } .hyp { background: var(--panel); padding: 10px 14px; border-radius: 8px; }
table { border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; } th, td { padding: 5px 8px; border-bottom: 1px solid var(--rule); text-align: left; }
th { font-weight: 600; color: var(--quiet); font-size: 12px; } td.better { color: var(--good); font-weight: 600; } td.worse { color: var(--bad); font-weight: 600; }
.key span { display: inline-block; width: 18px; height: 3px; vertical-align: middle; margin: 0 4px 0 12px; }
.snaps { display: grid; grid-template-columns: repeat(7, 1fr); gap: 6px; } .snap { margin: 0; } .snap img { width: 100%; border: 2px solid transparent; border-radius: 4px; background: #fff; }
.snap.breached img { border-color: var(--bad); } .snap figcaption { font-size: 11px; color: var(--quiet); text-align: center; }
.missing { aspect-ratio: 1; display: grid; place-items: center; color: var(--quiet); background: var(--panel); border-radius: 4px; }
.rowlabel { font-weight: 600; margin: 10px 0 4px; } .charts { display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 12px; }
.chart { width: 100%; height: auto; background: var(--panel); border-radius: 8px; } .ctitle { font-size: 12px; font-weight: 600; fill: var(--ink); }
.tick { font-size: 10px; fill: var(--quiet); } .grid { stroke: var(--rule); } .line-b { fill: none; stroke: var(--base); stroke-width: 1.8; }
.line-c { fill: none; stroke: var(--cand); stroke-width: 1.8; } .breach-b, .breach-c { stroke-dasharray: 4 3; stroke-width: 1.2; }
.breach-b { stroke: var(--base); } .breach-c { stroke: var(--cand); }
@media (max-width: 700px) { .snaps { grid-template-columns: repeat(4, 1fr); } }
"""

# ── The page ──
function build_report(outdir)
    sims = read_table(joinpath(outdir, "simulations.csv"))
    summary = read_table(joinpath(outdir, "summary.csv"))
    candidate_src = read(joinpath(outdir, "candidate.jl"), String)
    hypothesis = join([lstrip(l[2:end]) for l in split(candidate_src, "\n") if startswith(l, "#")], " ")
    hypothesis = replace(hypothesis, r"^Hypothesis:\s*"i => "")
    name = basename(outdir)
    experiments = unique(column(sims, "experiment"))
    sim_of(e, v) = Int(only(r["sim_id"] for r in sims if r["experiment"] == e && r["version"] == v))
    summ(e, v) = only(r for r in summary if r["experiment"] == e && r["version"] == v)
    rows_of(e, v) = read_table(joinpath(outdir, "checks", "$(e)_$(v).csv"))
    baseline_inputs = only(unique(r["inputs"] for r in sims if r["version"] == "baseline"))
    candidate_inputs = only(unique(r["inputs"] for r in sims if r["version"] == "candidate"))

    io = IOBuffer()
    print(io, "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
          "<title>Comparison: $(html(name))</title><style>$CSS</style></head><body><main>")

    # 1. Header
    print(io, "<h1>Comparison: $(html(name))</h1><p class=\"quiet\">Baseline <b>$(html(baseline_inputs))</b> vs candidate inputs <b>$(html(candidate_inputs))</b> · built $(Dates.format(now(), "yyyy-mm-dd HH:MM"))</p>",
          "<p class=\"hyp\"><b>Hypothesis.</b> $(html(hypothesis))</p>")
    print(io, "<h3>Settings that differ (first experiment, $(html(experiments[1])))</h3><table><tr><th>Setting</th><th>Baseline</th><th>Candidate</th></tr>")
    for (k, a, b) in settings_diff(sim_of(experiments[1], "baseline"), sim_of(experiments[1], "candidate"))
        print(io, "<tr><td>$(html(k))</td><td>$(html(a))</td><td>$(html(b))</td></tr>")
    end
    print(io, "</table>")

    # 2. Verdict table
    cols = [("fold_first", "Fold first?"), ("first_fold", "First fold (day)"), ("first_second_layer", "First 2nd layer (day)"),
            ("first_breach", "First breach (day)"), ("max_stretch", "Max stretch"), ("max_edge_stretch", "Max edge stretch"),
            ("min_curv_radius", "Min curvature radius (µm)"), ("max_disp", "Max displacement (µm)")]
    print(io, "<h2>Summary</h2><p class=\"quiet\">Stretch, curvature and displacement are worst values before the first breach. Green / red: candidate better / worse than baseline.</p>",
          "<table><tr><th>Experiment</th><th>Version</th>", join("<th>$h</th>" for (_, h) in cols), "</tr>")
    for e in experiments
        b, c = summ(e, "baseline"), summ(e, "candidate")
        print(io, "<tr><td rowspan=\"2\"><b>$(html(e))</b></td><td>baseline</td>", join("<td>$(fmt(b[k]))</td>" for (k, _) in cols), "</tr>")
        print(io, "<tr><td>candidate</td>", join("<td class=\"$(verdict(k, b[k], c[k]))\">$(fmt(c[k]))</td>" for (k, _) in cols), "</tr>")
    end
    print(io, "</table>")

    # Regression: identical to baseline?
    regfile = joinpath(outdir, "regression.csv")
    if isfile(regfile)
        print(io, "<h3>Identical to baseline?</h3><p class=\"quiet\">Run because the candidate sets <code>check_identical = true</code>; ",
              "a code cleanup or the negative control should read <b>yes</b>.</p><table><tr><th>Experiment</th><th>Identical</th><th>First difference</th></tr>")
        for r in read_table(regfile)
            print(io, "<tr><td>$(html(r["experiment"]))</td><td><b>$(fmt(r["identical"]))</b></td><td>$(fmt(r["first_difference"]))</td></tr>")
        end
        print(io, "</table>")
    end

    # 3 + 4. Per experiment: snapshots, then charts
    for e in experiments
        b_rows, c_rows = rows_of(e, "baseline"), rows_of(e, "candidate")
        b, c = summ(e, "baseline"), summ(e, "candidate")
        print(io, "<h2>$(html(e))</h2>")
        print(io, "<div class=\"rowlabel\">Baseline (sim $(sim_of(e, "baseline")))</div><div class=\"snaps\">", snapshot_row(sim_of(e, "baseline"), b_rows, b), "</div>")
        print(io, "<div class=\"rowlabel\">Candidate (sim $(sim_of(e, "candidate")))</div><div class=\"snaps\">", snapshot_row(sim_of(e, "candidate"), c_rows, c), "</div>")
        print(io, "<p class=\"key quiet\"><span style=\"background:var(--base)\"></span>baseline<span style=\"background:var(--cand)\"></span>candidate · dashed line = first breach</p><div class=\"charts\">")
        for (title, unit, col, zero) in CHARTS
            print(io, chart(title, unit, b_rows, c_rows, col, b["first_breach"], c["first_breach"]; from_zero = zero))
        end
        print(io, "</div>")
    end

    # 5. Notes
    print(io, "<h2>What the numbers mean</h2><table>", join("<tr><td><b>$(html(k))</b></td><td>$(html(v))</td></tr>" for (k, v) in NOTES), "</table></main></body></html>")

    path = joinpath(outdir, "report.html")
    write(path, take!(io))
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__   # standalone: rebuild from saved files
    using PhysiCellModelManager
    initializeModelManager(joinpath(@__DIR__, "..", ".."))
    include(joinpath(@__DIR__, "checks.jl"))   # output_dir
    isempty(ARGS) && error("name a candidate, e.g.  julia --project=. Model_Testing/report.jl adhesion_0.1")
    println(build_report(joinpath(@__DIR__, "Outputs", ARGS[1])))
end
