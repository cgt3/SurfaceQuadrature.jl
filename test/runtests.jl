using SurfaceQuadrature
using GaussQuadrature
using Test

@testset "SurfaceQuadrature.jl" begin
    @testset "get_line_quadrature" begin
        testing_tol = 1e-15

        # 1. Identity mapping
        @testset "Identity mapping" begin
            ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
            ref_pts_2D = [ [ref_pts[i], 0] for i in 1:length(ref_pts)]

            line(s) = [s, 0*s] # the x-axis
            s_domain = (-1, 1)

            line_pts, line_wts, normals = get_line_quadrature((ref_pts, ref_wts), line, s_domain)
            @test line_wts == ref_wts
            @test length(line_pts) == length(ref_pts)

            for i = 1:length(ref_pts)
                diff = abs.(line_pts[i] - ref_pts_2D[i])
                @test maximum(diff) < testing_tol
                @test normals[i] == [0, -1]
            end
        end

        # 2. Linear mapping to test normalization
        @testset "Linear mapping w normalization" begin
            ref_pts, ref_wts = legendre(5) # Legendre-Gauss quadrature
            true_pts_2D = [ [ref_pts[i], ref_pts[i]] for i in 1:length(ref_pts)]

            line(s) = [s, s] # simple linear function
            s_domain = (-1, 1)

            line_pts, line_wts, normals = get_line_quadrature((ref_pts, ref_wts), line, s_domain, normalization=true)
            @test line_wts == sqrt(2)*ref_wts
            @test length(line_pts) == length(ref_pts)

            for i = 1:length(ref_pts)
                diff_pts = abs.(line_pts[i] - true_pts_2D[i])
                @test maximum(diff_pts) < testing_tol
                
                diff_normals = abs.(normals[i] - sqrt(2)/2*[1, -1])
                @test maximum(diff_normals) < testing_tol
            end
        end


    end # test_set: "get_line_quadrature"
end
