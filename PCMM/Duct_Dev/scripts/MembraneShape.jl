#=
MembraneShape.jl — pure geometric shape descriptors for a closed 2D contour
============================================================================
This module is the formalized "statistics" core of the calibration pipeline.

Design principles:
  • PURE GEOMETRY, NO SIMULATION I/O. Every function takes an N×2 matrix of
    (x, y) vertices (a closed polygon, in order) and returns numbers. This keeps
    the descriptors unit-testable against analytic shapes (see TestDescriptors.jl)
    and independent of PhysiCell / PCMM.
  • H&E FEASIBILITY. Every descriptor is computable from a *single* traced
    contour using intrinsic (per-contour) references — no per-node identity and
    no initial-frame reference — so the same code can later run on a lumen
    boundary traced from an H&E image, not just simulation output.
  • RADIAL PROFILE r(θ). We represent the contour by its radius as a function of
    angle about the centroid, resampled onto a uniform θ grid. This is robust to
    remeshing (node count changing between frames) because it describes the
    *contour*, not indexed nodes.

Scale normalizer: R_eff = sqrt(area / π), the radius of the equal-area circle.

All descriptors returned by `shape_descriptors` are dimensionless (~O(0.1–1)),
so a mean-squared-error distance over them is well-behaved without ad-hoc scaling.
=#

using Statistics

# ──────────────────────────────────────────────────────────────────
#  BASIC POLYGON QUANTITIES  (operate on an N×2 matrix of vertices)
# ──────────────────────────────────────────────────────────────────

"""Signed polygon area via the shoelace formula (positive for counter-clockwise)."""
function polygon_signed_area(pts::AbstractMatrix{<:Real})::Float64
    N = size(pts, 1)
    a = 0.0
    @inbounds for i in 1:N
        j = (i % N) + 1
        a += pts[i, 1] * pts[j, 2] - pts[j, 1] * pts[i, 2]
    end
    return a / 2.0
end

"""Enclosed area of the polygon."""
polygon_area(pts::AbstractMatrix{<:Real})::Float64 = abs(polygon_signed_area(pts))

"""Perimeter — sum of consecutive edge lengths around the closed polygon."""
function polygon_perimeter(pts::AbstractMatrix{<:Real})::Float64
    N = size(pts, 1)
    p = 0.0
    @inbounds for i in 1:N
        j = (i % N) + 1
        p += hypot(pts[j, 1] - pts[i, 1], pts[j, 2] - pts[i, 2])
    end
    return p
end

"""
Area-weighted polygon centroid (robust to non-uniform vertex spacing, unlike a
plain vertex mean). Falls back to the vertex mean for a degenerate (~zero-area)
polygon.
"""
function polygon_centroid(pts::AbstractMatrix{<:Real})::Tuple{Float64,Float64}
    N = size(pts, 1)
    A = polygon_signed_area(pts)
    if abs(A) < 1e-12
        return (mean(@view pts[:, 1]), mean(@view pts[:, 2]))
    end
    cx = 0.0
    cy = 0.0
    @inbounds for i in 1:N
        j = (i % N) + 1
        cross = pts[i, 1] * pts[j, 2] - pts[j, 1] * pts[i, 2]
        cx += (pts[i, 1] + pts[j, 1]) * cross
        cy += (pts[i, 2] + pts[j, 2]) * cross
    end
    return (cx / (6A), cy / (6A))
end

"""Effective radius of a region with the given area: R_eff = sqrt(area / π)."""
effective_radius(area::Real)::Float64 = sqrt(area / π)

# ──────────────────────────────────────────────────────────────────
#  RADIAL PROFILE  r(θ)
# ──────────────────────────────────────────────────────────────────

"""
Periodic linear interpolation of samples (θs, rs) — θs sorted ascending in
[-π, π) — evaluated at query angle `tq` ∈ [-π, π). Wraps around the seam so the
first and last samples connect through ±π.
"""
function _interp_periodic(θs::AbstractVector{<:Real}, rs::AbstractVector{<:Real}, tq::Real)::Float64
    n = length(θs)
    if tq < θs[1]                       # seam: between θs[n]-2π and θs[1]
        θ0 = θs[n] - 2π; θ1 = θs[1]; r0 = rs[n]; r1 = rs[1]
    elseif tq >= θs[n]                  # seam: between θs[n] and θs[1]+2π
        θ0 = θs[n]; θ1 = θs[1] + 2π; r0 = rs[n]; r1 = rs[1]
    else
        idx = searchsortedfirst(θs, tq) # first θs[idx] >= tq, with idx >= 2 here
        θs[idx] == tq && return rs[idx]
        θ0 = θs[idx-1]; θ1 = θs[idx]; r0 = rs[idx-1]; r1 = rs[idx]
    end
    Δ = θ1 - θ0
    Δ == 0 && return r0
    return r0 + (tq - θ0) / Δ * (r1 - r0)
end

"""
    radial_profile(pts; M=256) -> (θgrid, r)

Resample the contour onto a uniform grid of `M` angles about its centroid.
Returns `θgrid` (length M, in [-π, π)) and `r` (radius at each θ). Handles
arbitrary vertex count/order and remeshing by sorting on angle and interpolating.
Assumes the contour is star-shaped about its centroid (true for the duct lumen
geometries here); multi-valued r(θ) from severe self-folding would be flattened.
"""
function radial_profile(pts::AbstractMatrix{<:Real}; M::Int=256)
    cx, cy = polygon_centroid(pts)
    N = size(pts, 1)
    θ = Vector{Float64}(undef, N)
    r = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        dx = pts[i, 1] - cx
        dy = pts[i, 2] - cy
        θ[i] = atan(dy, dx)
        r[i] = hypot(dx, dy)
    end
    perm = sortperm(θ)
    θs = θ[perm]
    rs = r[perm]
    θgrid = collect(range(-π, stop = π, length = M + 1))[1:M]
    rg = similar(θgrid)
    @inbounds for k in 1:M
        rg[k] = _interp_periodic(θs, rs, θgrid[k])
    end
    return θgrid, rg
end

# ──────────────────────────────────────────────────────────────────
#  FOURIER AMPLITUDES  (dependency-free; the pipeline avoids FFTW)
# ──────────────────────────────────────────────────────────────────

"""
Magnitude spectrum of a uniformly-sampled periodic signal `f` (length M).
Returns `amp` where `amp[k+1]` is the amplitude at integer wavenumber
k = 0 … M÷2. Computed by direct DFT sums (M small here, ~256), so no FFTW.
"""
function fourier_amplitudes(f::AbstractVector{<:Real})::Vector{Float64}
    M = length(f)
    nyq = M ÷ 2
    amp = zeros(Float64, nyq + 1)
    @inbounds for k in 0:nyq
        cr = 0.0
        ci = 0.0
        ω = 2π * k / M
        for n in 0:(M-1)
            cr += f[n+1] * cos(ω * n)
            ci -= f[n+1] * sin(ω * n)
        end
        amp[k+1] = hypot(cr, ci) * (2.0 / M)
    end
    return amp
end

# ──────────────────────────────────────────────────────────────────
#  DESCRIPTORS
# ──────────────────────────────────────────────────────────────────

"""
    shape_descriptors(pts; M=256, δ=0.05, lobe_band=2:6, rough_kmin=10)

Compute the normalized shape descriptors for one closed contour. All fields are
dimensionless and computable from a single traced contour (H&E-transferable),
except `area`/`perimeter`/`R_eff`, which are returned raw for reference.

Fields:
  IC                inverse circularity  P²/(4πA)                (1.0 = circle)
  max_indent_depth  maxθ (R_eff − r(θ)) / R_eff                  deepest inward dip
  indent_extent     fraction of θ-grid with depth(θ) > δ         width of deformation
  lobe_amp          max amplitude of r(θ)/R_eff over low wavenumbers `lobe_band`
                                                                  (focal / multi-lobe deformation)
  roughness         RMS amplitude of r(θ)/R_eff above wavenumber `rough_kmin`
                                                                  (jaggedness / instability signature)
  R_eff, area, perimeter                                         raw references
"""
function shape_descriptors(pts::AbstractMatrix{<:Real}; M::Int = 256, δ::Float64 = 0.05,
                           lobe_band::UnitRange{Int} = 2:6, rough_kmin::Int = 10)
    area = polygon_area(pts)
    perim = polygon_perimeter(pts)
    Reff = effective_radius(area)
    IC = area > 0 ? perim^2 / (4π * area) : NaN

    _, r = radial_profile(pts; M = M)
    rn = r ./ Reff                       # normalized radius  r(θ)/R_eff
    depth = 1.0 .- rn                    # (R_eff − r)/R_eff; positive = inward dip
    max_indent_depth = maximum(depth)
    indent_extent = count(>(δ), depth) / M

    # Spectrum of the normalized-radius deviation (DC removed).
    amp = fourier_amplitudes(rn .- mean(rn))
    nyq = length(amp) - 1
    amp_k(k) = (0 <= k <= nyq) ? amp[k+1] : 0.0

    klo = first(lobe_band)
    khi = min(last(lobe_band), nyq)
    lobe_amp = khi >= klo ? maximum(amp_k(k) for k in klo:khi) : 0.0

    # Roughness = RMS amplitude of the high-wavenumber (k > rough_kmin) radial
    # deviation, normalized by R_eff. Absolute (not a fraction of total energy),
    # so a smooth contour → ~0 rather than amplifying numerical noise into a
    # spurious large fraction.
    hi_energy = sum(amp_k(k)^2 for k in (rough_kmin+1):nyq; init = 0.0)  # k > rough_kmin
    roughness = sqrt(hi_energy)

    return (IC = IC,
            max_indent_depth = max_indent_depth,
            indent_extent = indent_extent,
            lobe_amp = lobe_amp,
            roughness = roughness,
            R_eff = Reff,
            area = area,
            perimeter = perim)
end
