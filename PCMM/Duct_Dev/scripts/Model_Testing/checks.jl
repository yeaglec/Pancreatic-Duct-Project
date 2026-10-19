# Checks computed from a simulation's saved outputs only (boundary_t<frame>.csv + cell positions), never model internals.
#
# Per saved frame (frame_checks):
#   Health      outside          cells whose centre is outside the membrane (escaped)
#               wall_gap         closest approach of two membrane points that are >= ARC_EXCLUDE apart along the
#                                membrane, i.e. two walls of a fold touching (µm)
#               longest_edge     longest membrane segment (µm)
#               disp_max/mean    distance of the membrane from its starting contour (µm); row 1 at rest
#   Numerical   stretch          perimeter / starting perimeter (BM should fold by bending, not stretching)
#               edge_stretch     largest edge length / its starting length; only while the node count is unchanged
#               min_curv_radius  smallest radius of a circle through three membrane points CURV_SPAN apart (µm);
#                                spikes show as radii far below a cell size, node jitter does not
#   Metrics.jl  IC, indent, off_membrane (fraction), filling, plus cell and off-membrane counts
#   Order       disp_since_division  max distance of the membrane from where it was at the first division (µm;
#                                    NaN before it), so startup settling doesn't count as a lesion fold
#
# Per simulation (summarize): first division, first fold, first second-layer cell, order of events, first breach
# (first frame of >= BREACH_FRAMES consecutive frames with a cell outside), and the worst stretch / curvature
# before the first breach. Startup effects are excluded from the order of events: the fold and the second layer
# are both measured from the first division.
# Between simulations (same_output): the regression check, i.e. are two runs identical?
using DelimitedFiles
include(joinpath(@__DIR__, "..", "Metrics.jl"))   # area, perimeter, distance_to_membrane, the four metrics

const R_CELL = 8.4          # cell radius (µm): the fold threshold and the curvature span
const ARC_EXCLUDE = 40.0    # µm along the membrane: closer pairs are neighbours, not opposite walls
const CURV_SPAN = R_CELL    # µm along the membrane between the three points of the curvature circle
const REACH = 10.0          # µm: membrane adhesion range, the off-membrane cutoff (as in the figures)
const BREACH_FRAMES = 2     # consecutive frames with a cell outside before it counts as a breach (1 = strict)

# ── Loading one saved frame ──
output_dir(sim_id) = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output")
load_membrane(sim_id, frame) = readdlm(joinpath(output_dir(sim_id), "boundary_t$(frame).csv"), ',', Float64)[:, 1:2]
function load_cells(sim_id, frame)
    c = PhysiCellSnapshot(sim_id, frame; include_cells = true).cells
    c = c[(.!c.dead) .& (c.cell_type_name .!= "BM"), :]
    return c.position_1, c.position_2, c.radius
end
n_frames(sim_id) = count(f -> startswith(f, "boundary_t") && endswith(f, ".csv"), readdir(output_dir(sim_id)))
function save_interval(sim_id)   # minutes between full-data saves, read from the settings the simulation ran with
    xml = read(joinpath(output_dir(sim_id), "PhysiCell_settings.xml"), String)
    return parse(Float64, match(r"<full_data>.*?<interval[^>]*>([^<]+)<"s, xml)[1])
end

# ── Geometry helpers ──
function inside(m, x, y)   # ray casting
    n, c = size(m, 1), false
    for i in 1:n
        j = i == 1 ? n : i - 1
        if (m[i, 2] > y) != (m[j, 2] > y) && x < (m[j, 1] - m[i, 1]) * (y - m[i, 2]) / (m[j, 2] - m[i, 2]) + m[i, 1]
            c = !c
        end
    end
    return c
end
edge_lengths(m) = [hypot(m[mod1(i + 1, size(m, 1)), 1] - m[i, 1], m[mod1(i + 1, size(m, 1)), 2] - m[i, 2]) for i in 1:size(m, 1)]
arclength(m) = cumsum([0.0; edge_lengths(m)])   # s[i] at node i; s[end] = perimeter

function wall_gap(m)
    n, s = size(m, 1), arclength(m)
    P, best = s[end], Inf
    for i in 1:n, j in (i + 1):n
        a = s[j] - s[i]
        min(a, P - a) < ARC_EXCLUDE && continue
        best = min(best, hypot(m[i, 1] - m[j, 1], m[i, 2] - m[j, 2]))
    end
    return best
end

function point_at(m, s, t)   # the membrane point at arc length t (wraps around)
    P = s[end]
    t = mod(t, P)
    k = clamp(searchsortedlast(s, t), 1, size(m, 1))
    j = mod1(k + 1, size(m, 1))
    f = (t - s[k]) / (s[k + 1] - s[k])
    return m[k, 1] + f * (m[j, 1] - m[k, 1]), m[k, 2] + f * (m[j, 2] - m[k, 2])
end

function min_curvature_radius(m)
    s, best = arclength(m), Inf
    for i in 1:size(m, 1)
        ax, ay = point_at(m, s, s[i] - CURV_SPAN)
        bx, by = m[i, 1], m[i, 2]
        cx, cy = point_at(m, s, s[i] + CURV_SPAN)
        a, b, c = hypot(bx - cx, by - cy), hypot(ax - cx, ay - cy), hypot(ax - bx, ay - by)
        twice_area = abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))
        twice_area > 1e-12 && (best = min(best, a * b * c / (2 * twice_area)))   # circumradius abc / 4A
    end
    return best
end

# ── Per frame ──
function frame_checks(sim_id, frame, m0)
    m = load_membrane(sim_id, frame)
    x, y, r = load_cells(sim_id, frame)
    disp = [distance_to_membrane(m0, m[i, 1], m[i, 2]) for i in 1:size(m, 1)]
    off = off_membrane(m, x, y, r; reach = REACH)
    e, e0 = edge_lengths(m), edge_lengths(m0)
    return (outside = count(k -> !inside(m, x[k], y[k]), eachindex(x)),
            wall_gap = wall_gap(m), longest_edge = maximum(e),
            disp_max = maximum(disp), disp_mean = sum(disp) / length(disp),
            stretch = perimeter(m) / perimeter(m0),
            edge_stretch = length(e) == length(e0) ? maximum(e ./ e0) : NaN,
            min_curv_radius = min_curvature_radius(m),
            IC = inverse_circularity(m), indent = max_indent_depth(m),
            off_membrane = sum(off) / length(x), filling = lumen_filling(m, r),
            cells = length(x), off_count = sum(off))
end

# All frames of one simulation as a vector of NamedTuples (one per frame), with time in days
function sim_checks(sim_id)
    m0, dt = load_membrane(sim_id, 0), save_interval(sim_id)
    rows = [merge((frame = f, day = f * dt / 1440), frame_checks(sim_id, f, m0)) for f in 0:(n_frames(sim_id) - 1)]
    k = findfirst(r -> r.cells > rows[1].cells, rows)   # first division
    k === nothing && return [merge(r, (disp_since_division = NaN,)) for r in rows]
    m_div = load_membrane(sim_id, rows[k].frame)
    since(r) = (m = load_membrane(sim_id, r.frame); maximum(distance_to_membrane(m_div, m[i, 1], m[i, 2]) for i in 1:size(m, 1)))
    return [merge(r, (disp_since_division = r.frame < rows[k].frame ? NaN : since(r),)) for r in rows]
end

# ── Per simulation ──
first_frame(rows, test) = (k = findfirst(test, rows); k === nothing ? missing : rows[k].day)

function summarize(rows)
    n0 = rows[1].cells
    k_div = findfirst(r -> r.cells > n0, rows)
    # second layer: off-membrane count rises above its value at the first division (ignores startup crowding)
    second = k_div === nothing ? missing : first_frame(rows[k_div:end], r -> r.off_count > rows[k_div].off_count)
    fold = first_frame(rows, r -> r.disp_since_division >= R_CELL)   # NaN before the division never passes
    n = BREACH_FRAMES - 1
    breach = findfirst(k -> k + n <= length(rows) && all(r -> r.outside > 0, rows[k:k + n]), eachindex(rows))
    valid = breach === nothing ? rows : rows[1:breach - 1]   # frames before the first escape
    return (first_division = k_div === nothing ? missing : rows[k_div].day,
            first_fold = fold, first_second_layer = second,
            fold_first = k_div === nothing ? missing : ismissing(second) ? true : (!ismissing(fold) && fold <= second),
            first_breach = breach === nothing ? missing : rows[breach].day,
            max_stretch = maximum(r -> r.stretch, valid), max_edge_stretch = maximum(r -> r.edge_stretch, valid),
            min_curv_radius = minimum(r -> r.min_curv_radius, valid), max_disp = maximum(r -> r.disp_max, valid))
end

# ── Regression: are two simulations identical? ──
# Boundary files byte for byte, cell positions as exact numbers (not the .mat files, which may carry metadata).
function same_output(sim_a, sim_b)
    n_frames(sim_a) == n_frames(sim_b) || return (identical = false, first_difference = "frame count")
    for f in 0:(n_frames(sim_a) - 1)
        read(joinpath(output_dir(sim_a), "boundary_t$f.csv")) == read(joinpath(output_dir(sim_b), "boundary_t$f.csv")) ||
            return (identical = false, first_difference = "membrane frame $f")
        load_cells(sim_a, f) == load_cells(sim_b, f) || return (identical = false, first_difference = "cells frame $f")
    end
    return (identical = true, first_difference = "")
end

# ── Writing ──
function write_csv(path, rows)
    open(path, "w") do io
        println(io, join(keys(first(rows)), ","))
        for r in rows
            println(io, join((v isa AbstractFloat ? round(v; sigdigits = 6) : v for v in values(r)), ","))
        end
    end
end
