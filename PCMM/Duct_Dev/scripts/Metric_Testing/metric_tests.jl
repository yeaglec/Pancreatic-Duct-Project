# Is each lesion metric correct? Synthetic membranes and cells with known answers; no simulations, runs in seconds.
# Usage (from scripts/):  julia --project=. Metric_Testing/metric_tests.jl
using Test
include(joinpath(@__DIR__, "..", "Metrics.jl"))

const D = D_CELL      # cell diameter
const R = D_CELL / 2  # cell radius

# ── Synthetic membranes ──
# Polygon through the given corners, with a point every ~1 µm
function polygon(corners)
    n = size(corners, 1)
    pts = [corners[i, :] + t * (corners[mod1(i + 1, n), :] - corners[i, :])
           for i in 1:n for t in range(0, 1; length = ceil(Int, hypot((corners[mod1(i + 1, n), :] - corners[i, :])...)) + 1)[1:end - 1]]
    return permutedims(reduce(hcat, pts))
end
# 300 µm square with a notch (depth d, width w) cut inward from the middle of the top edge (top at y = 150)
notched_square(d, w) = polygon([-150 -150; 150 -150; 150 150; w/2 150; w/2 150-d; -w/2 150-d; -w/2 150; -150 150])
# Circle of radius 150 with a smooth inward dent of depth d at the top, drawn with n nodes
dented_circle(d; n = 600) = permutedims(reduce(hcat, (150 - d * exp(-atan(sin(t - π/2), cos(t - π/2))^2 / (2 * 0.15^2))) * [cos(t), sin(t)]
                                                    for t in range(0, 2π; length = n + 1)[1:n]))
circle() = dented_circle(0.0)

# ── Synthetic cells ──
# k cells side by side at the top of the circle, centres `inset` µm inside the circle
function row(k, inset)
    θ = [π/2 + (j - (k + 1) / 2) * D / 150 for j in 1:k]
    return (150 - inset) .* cos.(θ), (150 - inset) .* sin.(θ)
end
metrics(m, x, y) = lesion_metrics(m, x, y, trues(length(x)))   # every cell is a lesion cell
rotate(m, a) = m * [cos(a) sin(a); -sin(a) cos(a)]

@testset "Lesion metrics" begin
    @testset "circle, one layer: nothing happened" begin
        g = metrics(circle(), row(5, R)...)
        @test g.IC ≈ 1 atol = 1e-3
        @test g.fold_depth ≈ 0 atol = 0.01
        @test g.fold_aspect == 0
        @test g.layer_thickness ≈ 1 atol = 0.02
        @test isnan(g.fold_share)
    end

    @testset "circle, two layers: stacking only" begin
        (x1, y1), (x2, y2) = row(5, R), row(5, R + D)
        g = metrics(circle(), [x1; x2], [y1; y2])
        @test g.layer_thickness ≈ 2 atol = 0.02
        @test g.fold_share ≈ 0 atol = 0.01
    end

    @testset "notch 2 deep, 2 wide: depth 2, aspect 1" begin
        g = metrics(notched_square(2D, 2D), [0.0], [150 - 2D - R])   # one cell at the bottom
        @test g.fold_depth ≈ 2 atol = 0.02
        @test g.fold_aspect ≈ 1 atol = 0.15                          # mouth resolution ≈ the 2 µm spacing
    end

    @testset "notch 3 deep, 1 wide, lined by cells: folding only" begin
        g = metrics(notched_square(3D, D), zeros(3), [150 - 3D + R + j * D for j in 0:2])
        @test g.fold_depth ≈ 3 atol = 0.02
        @test g.fold_aspect ≈ 3 rtol = 0.2
        @test g.layer_thickness ≈ 1 atol = 0.02                      # a cell deep in a fold is still on the membrane
        @test g.fold_share ≈ 1 atol = 0.01
    end

    @testset "dent shallower than MOUTH_EPS: aspect 0, not Inf" begin
        @test metrics(dented_circle(0.5), row(3, 0.5 + R)...).fold_aspect == 0
    end

    @testset "local: a fold away from the lesion is ignored" begin
        @test metrics(notched_square(2D, 2D), [0.0], [-150 + R]).fold_depth ≈ 0 atol = 0.01
    end

    @testset "rotation does not matter" begin
        m, cell = notched_square(2D, 2D), [0.0 150 - 2D - R]
        a = metrics(m, cell[:, 1], cell[:, 2])
        b = metrics(rotate(m, 0.5), rotate(cell, 0.5)[:, 1], rotate(cell, 0.5)[:, 2])
        @test all(isapprox(a[k], b[k]; atol = 0.02) for k in keys(a))
    end

    @testset "node spacing does not matter (0.5 vs 5 µm)" begin
        x, y = row(3, 25 + R)
        a, b = metrics(dented_circle(25; n = 1900), x, y), metrics(dented_circle(25; n = 190), x, y)
        @test all(isapprox(a[k], b[k]; rtol = 0.03) for k in keys(a))
    end
end
