#=
TestDescriptors.jl — analytic ground-truth tests for MembraneShape.jl
=====================================================================
These are the PRIMARY correctness check for the shape-descriptor module. They
run with no simulation and no PCMM — pure Julia + stdlib Test — by feeding
analytic contours (circle, ellipse, dented circle, cos(kθ)-perturbed circle)
whose descriptors are known in closed form.

Run from the scripts directory:
    julia TestDescriptors.jl
=#

using Test
include(joinpath(@__DIR__, "MembraneShape.jl"))

# ── analytic contour generators ────────────────────────────────────
circle(R; n = 400) = _polar(θ -> R, n)
ellipse(a, b; n = 400) = reduce(vcat, [a*cos(θ) b*sin(θ)] for θ in _angles(n))
cos_perturbed(R, ε, freq; n = 400) = _polar(θ -> R * (1 + ε * cos(freq * θ)), n)

function dented_circle(R, d0, θ0, width; n = 400)
    _polar(θ -> R - d0 * exp(-(_angdiff(θ, θ0)^2) / (2 * width^2)), n)
end

_angles(n) = [2π * i / n for i in 0:(n-1)]
_polar(rfun, n) = reduce(vcat, ([rfun(θ)*cos(θ) rfun(θ)*sin(θ)] for θ in _angles(n)))
_angdiff(a, b) = atan(sin(a - b), cos(a - b))

# ── tests ──────────────────────────────────────────────────────────
@testset "MembraneShape descriptors" begin

    @testset "perfect circle" begin
        d = shape_descriptors(circle(100.0))
        @test isapprox(d.IC, 1.0; atol = 1e-3)
        @test isapprox(d.R_eff, 100.0; atol = 0.5)
        @test isapprox(d.max_indent_depth, 0.0; atol = 1e-3)
        @test d.indent_extent == 0.0
        @test isapprox(d.lobe_amp, 0.0; atol = 1e-3)
        @test isapprox(d.roughness, 0.0; atol = 1e-3)
    end

    @testset "ellipse is a k=2 shape, IC>1" begin
        pts = ellipse(130.0, 100.0)
        d = shape_descriptors(pts)
        @test d.IC > 1.0
        @test d.lobe_amp > 0.0
        _, r = radial_profile(pts)
        amp = fourier_amplitudes(r ./ mean(r))
        @test argmax(amp[2:end]) == 2          # dominant wavenumber is k=2 (period π)
    end

    @testset "cos(kθ) perturbation → peak at k, amplitude ≈ ε" begin
        for (ε, freq) in [(0.08, 4), (0.05, 3), (0.10, 6)]
            pts = cos_perturbed(100.0, ε, freq)
            _, r = radial_profile(pts)
            amp = fourier_amplitudes(r ./ mean(r))
            @test argmax(amp[2:end]) == freq                 # dominant wavenumber is `freq`
            @test isapprox(amp[freq+1], ε; atol = 0.01)      # amplitude recovers ε
            d = shape_descriptors(pts)
            @test isapprox(d.lobe_amp, ε; atol = 0.015)      # (freq within default lobe band 2:6)
        end
    end

    @testset "single Gaussian dent → depth and extent" begin
        R, d0, width = 100.0, 20.0, 0.15
        d = shape_descriptors(dented_circle(R, d0, 0.7, width))
        # normalized depth uses R_eff (< R here, since the dent removes area) and
        # the narrow dent bottom is slightly under-resolved, so it lands a bit
        # below the naive d0/R = 0.2.
        @test 0.14 < d.max_indent_depth < 0.20
        @test 0.0 < d.indent_extent < 0.5                          # focal, not whole ring

        # wider dent ⇒ larger angular extent (monotonicity)
        narrow = shape_descriptors(dented_circle(R, d0, 0.7, 0.10)).indent_extent
        wide   = shape_descriptors(dented_circle(R, d0, 0.7, 0.30)).indent_extent
        @test wide > narrow

        # deeper dent ⇒ larger max_indent_depth
        shallow = shape_descriptors(dented_circle(R, 10.0, 0.7, width)).max_indent_depth
        deep    = shape_descriptors(dented_circle(R, 30.0, 0.7, width)).max_indent_depth
        @test deep > shallow
    end

    @testset "roughness detects high-frequency jaggedness" begin
        smooth = shape_descriptors(cos_perturbed(100.0, 0.08, 3)).roughness
        jagged = shape_descriptors(cos_perturbed(100.0, 0.08, 20)).roughness   # k=20 > rough_kmin
        @test jagged > smooth
        @test smooth < 0.1
    end

    @testset "invariance: translation & scale" begin
        base = cos_perturbed(100.0, 0.08, 4)
        shifted = base .+ [37.0 -22.0]                 # translate
        scaled = 2.5 .* base                           # scale
        db = shape_descriptors(base)
        ds = shape_descriptors(shifted)
        dc = shape_descriptors(scaled)
        for f in (:IC, :max_indent_depth, :indent_extent, :lobe_amp, :roughness)
            @test isapprox(getfield(db, f), getfield(ds, f); atol = 1e-6)   # translation-invariant
            @test isapprox(getfield(db, f), getfield(dc, f); atol = 1e-6)   # scale-invariant
        end
    end
end
