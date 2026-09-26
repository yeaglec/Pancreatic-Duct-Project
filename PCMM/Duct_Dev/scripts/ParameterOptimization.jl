#=
ParameterOptimization.jl — ABC-SMC Calibration for Pancreatic Duct Models
===========================================================================
Provides custom basement membrane (BM) morphology metrics and wraps them
into PCMM's calibration interface.

PCMM's calibration module handles all sampling, perturbation, ect. internally. We supply just two things:
  1. A summary statistic:  (monad_id) → Dict{String, Float64}
  2. A distance function:  (simulated, observed) → Float64

Upgrade PCMM to v0.3.0+ (from main) to use the calibration module:
    using Pkg
    Pkg.add(url="https://github.com/drbergman-lab/PhysiCellModelManager.jl", rev="main")

=#

using PhysiCellModelManager
using DelimitedFiles
using Statistics

# Pure-geometry shape descriptors (radial-profile based, H&E-transferable).
# This is the formalized statistics core; keep sim I/O here, geometry there.
include(joinpath(@__DIR__, "MembraneShape.jl"))

# ══════════════════════════════════════════════════════════════════
#                  BASEMENT MEMBRANE METRICS
# ══════════════════════════════════════════════════════════════════
# These operate on boundary_t{index}.csv files produced by the simulation.

"""Output path for a given simulation ID."""
_sim_output_dir(sim_id::Int) = joinpath(
    PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output")

"""Return sorted list of available boundary timestep indices."""
function boundary_timesteps(sim_id::Int)::Vector{Int}
    out_dir = _sim_output_dir(sim_id)
    !isdir(out_dir) && error("Output dir not found: $out_dir")

    # readdir - gets a list of every file in outdir
    # filter - loops through this list and keeps only boundary.csv files
    # Look at only boundary.csv files
    files = filter(f -> startswith(f, "boundary_t") && endswith(f, ".csv"), readdir(out_dir))       # -> syntax is for lambda function, which the filter only includes if evaluated to true

    # parse - converts string to type int
    # match - check string for match to a regex, then returns the regex match
    # Stores only the boundary files needed
    ts = [parse(Int, m.captures[1]) for f in files
          for m in (match(r"boundary_t(\d+)\.csv$", f),) if m !== nothing]

    return sort(ts)
end

"""Load boundary as N×2 matrix. Defaults to final timestep."""
function _load_boundary(sim_id::Int; timestep::Union{Nothing,Int}=nothing)::Matrix{Float64}
    # Default: Last simulation output
    if timestep === nothing
        ts = boundary_timesteps(sim_id)
        isempty(ts) && error("No boundary CSVs for simulation $sim_id")
        timestep = last(ts)
    end

    csv_path = joinpath(_sim_output_dir(sim_id), "boundary_t$(timestep).csv")
    !isfile(csv_path) && error("File not found: $csv_path")
    return readdlm(csv_path, ',', Float64)
end

"""Perimeter of the BM polygon (sum of consecutive edge lengths)."""
Perimeter(sim_id::Int; timestep::Union{Nothing,Int}=nothing)::Float64 =
    polygon_perimeter(_load_boundary(sim_id; timestep))

"""Enclosed area of the BM polygon (Shoelace formula)."""
Area(sim_id::Int; timestep::Union{Nothing,Int}=nothing)::Float64 =
    polygon_area(_load_boundary(sim_id; timestep))

"""Area-to-perimeter ratio. For a circle of radius r this equals r/2."""
function AreaToPerimeterRatio(sim_id::Int; timestep::Union{Nothing,Int}=nothing)::Float64
    a = Area(sim_id; timestep)
    p = Perimeter(sim_id; timestep)
    p == 0.0 && error("Perimeter is zero for simulation $sim_id")
    return a / p
end

"""Inverse circularity: 1.0 for a circle, >1 for deformed shapes."""
function InverseCircularity(sim_id::Int; timestep::Union{Nothing,Int}=nothing)::Float64
    a = Area(sim_id; timestep)
    p = Perimeter(sim_id; timestep)
    a == 0.0 && error("Area is zero for simulation $sim_id")
    return p^2 / (4π * a)
end

"""Max displacement of any BM node from its initial position."""
function MaxNodeDisplacement(sim_id::Int; timestep::Union{Nothing,Int}=nothing)::Float64
    ts = boundary_timesteps(sim_id)
    isempty(ts) && error("No boundary CSVs for simulation $sim_id")

    pts_initial = _load_boundary(sim_id; timestep=first(ts))
    pts_final   = _load_boundary(sim_id; timestep=timestep)

    # Handle node count mismatch (from add_membrane_nodes)
    N = min(size(pts_initial, 1), size(pts_final, 1))
    max_d = 0.0
    for i in 1:N
        dx = pts_final[i,1] - pts_initial[i,1]
        dy = pts_final[i,2] - pts_initial[i,2]
        max_d = max(max_d, sqrt(dx^2 + dy^2))
    end
    return max_d
end

# TODO: Add statistics involving agents and not just BM
# Add statistics that account for local geometry (max deformation depth & width, lobe amplitude, roughness, etc.)
# CellBreachCount — count EP cells outside the BM polygon
# Others

# ══════════════════════════════════════════════════════════════════
#             CELL-COUNT MILESTONE → OUTPUT FRAME
# ══════════════════════════════════════════════════════════════════
# We enable descriptors  to be evaluated at a fixed cell-count milestone rather than the final frame.
# Shared with GenerateReport.jl's :by_cell_count snapshot selection so both pick the same frame.

"""
    output_index_at_cell_count(sim_id, target) -> Union{Nothing, Tuple{Int,Int}}

Return the first output snapshot where total cell count reaches `target`, as
`(index, count)` with a **0-based** index to match PCMM's snapshot/boundary
numbering. Return `nothing` if `target` is never reached. """
function output_index_at_cell_count(sim_id::Int, target::Int)
    spts = PhysiCellModelManager.SimulationPopulationTimeSeries(sim_id)
    for i in 1:length(spts.time)
        total = sum(counts[i] for counts in values(spts.cell_count); init=0)
        if total >= target
            return (i - 1, total)   # snapshot/boundary indices are 0-based
        end
    end
    return nothing
end

"""
    boundary_timestep_at_cell_count(sim_id, target) → Int

Boundary CSV timestep index for the cell-count milestone `target`. Falls back to
the final available boundary timestep if the milestone is never reached, the
series is unavailable, or the mapped index has no boundary file. """
function boundary_timestep_at_cell_count(sim_id::Int, target::Int)::Int
    ts = boundary_timesteps(sim_id)
    isempty(ts) && error("No boundary CSVs for simulation $sim_id")
    res = try
        output_index_at_cell_count(sim_id, target)
    catch e
        @warn "Could not load population time series for sim $sim_id; using final frame" exception=e
        nothing
    end
    if res === nothing
        # Milestone unreachable → the milestone eval silently becomes a final-frame
        # eval, so warn loudly (set target_cell_count to a value the sims reach).
        @warn "Sim $sim_id never reached target cell count $target; using final boundary frame $(last(ts))"
        return last(ts)
    end
    idx = first(res)
    idx in ts && return idx
    @warn "Mapped output index $idx has no boundary file for sim $sim_id; using final frame $(last(ts))"
    return last(ts)
end

# ══════════════════════════════════════════════════════════════════
#             PER-SIMULATION METRIC EVALUATION
# ══════════════════════════════════════════════════════════════════

"""
    evaluate_simulation(sim_id; target_cell_count=nothing, timestep=nothing) → Dict{String,Float64}

Compute all BM metrics for one simulation. The evaluation frame `tf` is resolved in priority order: an explicit `timestep`, else the `target_cell_count` milestone
(via `boundary_timestep_at_cell_count`), else the final frame.

Emits both the original keys (`IC`, `area_frac_change`, `max_displacement`, …) and the normalized shape descriptors from `MembraneShape.jl`
(`max_indent_depth`, `indent_extent`, `lobe_amp`, `roughness`). """
function evaluate_simulation(sim_id::Int; target_cell_count::Union{Nothing,Int}=nothing,
                             timestep::Union{Nothing,Int}=nothing)::Dict{String,Float64}
    # TODO: Support trajectory-based evaluation (multiple timesteps / steady state)

    ts = boundary_timesteps(sim_id)
    isempty(ts) && return Dict{String,Float64}()

    t0 = first(ts)
    tf = if timestep !== nothing
        timestep                                       # explicit frame wins
    elseif target_cell_count !== nothing
        boundary_timestep_at_cell_count(sim_id, target_cell_count)
    else
        last(ts)                                       # final frame by default
    end

    # Dictionary with metrics for the first (t0) and evaluation (tf) frame
    metrics = Dict{String,Float64}()
    try
        metrics["perimeter_t0"]  = Perimeter(sim_id; timestep=t0)
        metrics["area_t0"]       = Area(sim_id; timestep=t0)
        metrics["IC_t0"]         = InverseCircularity(sim_id; timestep=t0)
        metrics["perimeter_tf"]  = Perimeter(sim_id; timestep=tf)
        metrics["area_tf"]       = Area(sim_id; timestep=tf)
        metrics["IC_tf"]         = InverseCircularity(sim_id; timestep=tf)
        metrics["max_displacement"] = MaxNodeDisplacement(sim_id; timestep=tf)
        metrics["area_frac_change"] = (metrics["area_tf"] - metrics["area_t0"]) / metrics["area_t0"]
        metrics["eval_timestep"] = Float64(tf)

        # Normalized shape descriptors at the evaluation frame (single load).
        # These are the H&E-transferable, complementary local/global metrics.
        desc = shape_descriptors(_load_boundary(sim_id; timestep=tf))
        metrics["max_indent_depth"] = desc.max_indent_depth
        metrics["indent_extent"]    = desc.indent_extent
        metrics["lobe_amp"]         = desc.lobe_amp
        metrics["roughness"]        = desc.roughness

        # Summary keys expected by bm_summary_statistic / bm_distance
        metrics["IC"] = metrics["IC_tf"]
    catch e
        @warn "Metric evaluation failed for sim $sim_id" exception=e
    end

    return metrics
end

"""
    boundary_trajectory(sim_id; frames = :all) -> (times, series)

Shape descriptors across the boundary frames of ONE simulation.
  times  :: Vector{Int}          the boundary indices actually used
  series :: Vector{NamedTuple}   shape_descriptors(...) at each of those frames
"""
function boundary_trajectory(sim_id::Int; frames = :all)
    ts = boundary_timesteps(sim_id)
    isempty(ts) && return (Int[], NamedTuple[])

    sel = frames === :all ? ts : intersect(ts, frames)

    # Add additional metrics from shape_descriptors to the trajectory series
    series = map(sel) do t
        pts = _load_boundary(sim_id; timestep = t)   # load the frame ONCE
        A = polygon_area(pts)
        P = polygon_perimeter(pts)
        (IC = P^2 / (4π * A), A2P = A / P)           # the two metrics you're testing with

        # desc = shape_descriptors(pts)
    end
    return (sel, series)
end

"""
    tail_drift(values; window = 20) -> Float64

Mean absolute frame-to-frame change of a scalar time series over the last `window` frames. Returns `Inf` for series of length < 2.
"""
function tail_drift(values::AbstractVector{<:Real}; window::Int = 20)::Float64
    n = length(values)
    n < 2 && return Inf
    w = min(window, n)                       # clamp: use what we have if the sim is short
    tail = @view values[(n - w + 1):n]
    return mean(abs(tail[i+1] - tail[i]) for i in 1:(length(tail) - 1))
end

"""
    stability_metrics(sim_id; window = 20) -> Dict{String,Float64}

Trajectory-based stability summary for ONE simulation. Reduces the IC
time-course over the last `window` boundary frames to:

  terminal_drift  mean |ΔIC| over the window        ≈0 ⇒ shape has stopped moving (settled)
  IC              mean IC over the window            the shape it settled AT (anti-freezing target)
  n_frames        number of boundary frames          diagnostic (did the sim run to completion?)
"""
function stability_metrics(sim_id::Int; window::Int = 20)::Dict{String,Float64}
    _, series = boundary_trajectory(sim_id)
    isempty(series) && return Dict{String,Float64}()

    IC_t = getfield.(series, :IC)            # IC at every frame

    # WINDOW = fixed count of tail frames (~last 20 of ~96–192). Cheap: no
    # per-frame population lookup.
    # TODO Replicate :by_cell_count from GenerateReport (score stability only on
    # frames AFTER a cell-count threshold, i.e. under high load) instead of a
    # fixed tail. Needs SimulationPopulationTimeSeries per sim, which profiled as
    # too expensive across ~1000 sims — cache/precompute the cell-count series
    # (e.g. once per monad) before switching.
    w    = min(window, length(IC_t))         # clamp for short/crashed sims
    tail = @view IC_t[(end - w + 1):end]

    metrics = Dict{String,Float64}()
    metrics["terminal_drift"] = tail_drift(IC_t; window = window)  # ≈0 ⇒ settled
    metrics["IC"]             = mean(tail)                         # shape it settled AT
    metrics["n_frames"]       = Float64(length(series))
    return metrics
end

# ══════════════════════════════════════════════════════════════════
#             PCMM ABC-SMC CALIBRATION INTERFACE
# ══════════════════════════════════════════════════════════════════
# These two functions plug directly into CalibrationProblem.
#
# CalibrationProblem(ref, parameters, observed_data,
#     bm_summary_statistic,   # ← our summary statistic
#     bm_distance;            # ← our distance function
#     n_replicates = 2,
# )

"""
    bm_summary_statistic(monad_id::Int) → Dict{String, Float64}

Custom summary statistic for PCMM's ABC-SMC calibration. Called once per proposed particle. Retrieves all simulation IDs
from the monad, evaluates BM metrics on each, and averages across replicates. """
function bm_summary_statistic(monad_id::Int; target_cell_count::Union{Nothing,Int}=nothing)::Dict{String,Float64}
    sim_ids = simulationIDs(Monad(monad_id))

    if isempty(sim_ids)
        @warn "No simulations found for monad $monad_id"
        return Dict{String,Float64}()
    end

    # Evaluate each replicate (at the cell-count milestone when given)
    all_metrics = [evaluate_simulation(sid; target_cell_count=target_cell_count) for sid in sim_ids]

    # Filter out any that failed (empty dicts)
    valid = filter(!isempty, all_metrics)
    isempty(valid) && return Dict{String,Float64}()

    # Average across replicates
    all_keys = union([keys(m) for m in valid]...)
    avg = Dict{String,Float64}()
    for k in all_keys
        vals = [m[k] for m in valid if haskey(m, k)]
        avg[k] = isempty(vals) ? 0.0 : mean(vals)
    end

    return avg
end

"""
    stability_summary_statistic(monad_id; window = 20) → Dict{String, Float64}

Monad-level summary for the equilibrium/stability calibration. Mirrors
`bm_summary_statistic`, but reduces each replicate with the trajectory-based
`stability_metrics` (terminal_drift + settled IC) instead of the single-frame
`evaluate_simulation`, then averages across the monad's replicate sims.

NOTE on averaging: if any replicate blew up / died early, its `terminal_drift`
is `Inf`, so the averaged `terminal_drift` becomes `Inf` and the whole monad is
rejected by `bm_distance`. That is intentional — a parameter set that sometimes
explodes is not robustly stable. If that is too strict for seed noise, swap
`mean` for `median` below or drop the `Inf`s before averaging. """
function stability_summary_statistic(monad_id::Int; window::Int = 20)::Dict{String,Float64}
    sim_ids = simulationIDs(Monad(monad_id))

    if isempty(sim_ids)
        @warn "No simulations found for monad $monad_id"
        return Dict{String,Float64}()
    end

    # Reduce each replicate to its trajectory-based stability summary
    all_metrics = [stability_metrics(sid; window=window) for sid in sim_ids]

    valid = filter(!isempty, all_metrics)
    isempty(valid) && return Dict{String,Float64}()

    # Average across replicates (Inf in any replicate ⇒ Inf ⇒ monad rejected)
    all_keys = union([keys(m) for m in valid]...)
    avg = Dict{String,Float64}()
    for k in all_keys
        vals = [m[k] for m in valid if haskey(m, k)]
        avg[k] = isempty(vals) ? 0.0 : mean(vals)
    end

    return avg
end

# Per-metric scale s_k standardizes each residual so metrics on different scales
# contribute comparably (else IC ~O(1) swamps area_frac_change ~O(0.1)); weight
# w_k is optional emphasis on top. Missing keys fall back to 1.0.
# NEXT STEP: estimate scales from the LHS pilot's per-metric std (data-driven).
#
# terminal_drift: residuals are tiny (IC change per frame ~1e-3..1e-2) vs IC
# residuals (~0.1..0.3); without a matched scale the IC term swamps stability and
# the fit ignores drift. The 0.005 below is a PLACEHOLDER — set it to the
# steady-state drift noise floor (mean |ΔIC| measured on sims you know are stable).
const BM_SCALES  = Dict("IC" => 0.4, "area_frac_change" => 0.1, "terminal_drift" => 0.005)
const BM_WEIGHTS = Dict{String,Float64}()   # empty ⇒ equal (1.0) weight

"""
    bm_distance(simulated, observed) → Float64

Scale-normalized weighted MSE between simulated and observed BM metrics. Named
(not a closure) so the CalibrationProblem stays serializable for `resumeABC`."""
function bm_distance(simulated::Dict{String,Float64}, observed::Dict{String,Float64})::Float64
    return bm_weighted_distance(simulated, observed, BM_WEIGHTS, BM_SCALES)
end

"""
    d = (1/n) Σ_k  w_k · ((sim_k − obs_k) / s_k)²   over the keys in `observed`.

A missing key (failed sim) returns `Inf` — a clean ABC rejection."""
function bm_weighted_distance(simulated::Dict{String,Float64},
                              observed::Dict{String,Float64},
                              weights::Dict{String,Float64} = Dict{String,Float64}(),
                              scales::Dict{String,Float64}  = Dict{String,Float64}())::Float64
    d = 0.0
    n = 0
    for (key, obs_val) in observed
        haskey(simulated, key) || return Inf     # failed sim ⇒ reject
        s = get(scales, key, 1.0)
        w = get(weights, key, 1.0)
        d += w * ((simulated[key] - obs_val) / s)^2
        n += 1
    end
    return n > 0 ? d / n : Inf
end
