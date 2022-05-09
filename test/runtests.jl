using SurfaceQuadrature
using GaussQuadrature
using Test

@testset "SurfaceQuadrature.jl" begin
    @testset "get_line_quadrature" begin
        # 1. Identity mapping
        ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
        ref_pts_2D = [ [ref_pts[i], 0] for i in 1:length(ref_pts)]

        line(s) = [s, 0*s] # simple linear function
        s_domain = (-1, 1)

        line_pts, line_wts, normals = get_line_quadrature((ref_pts, ref_wts), line, s_domain)
        @test line_wts == ref_wts
        @test length(line_pts) == length(ref_pts)

        for i = 1:length(ref_pts)
            diff = abs.(line_pts[i] - ref_pts_2D[i])
            @test maximum(diff) < 1e-15
        end

    end # "get_line_quadrature"
end
