module SurfaceQuadrature

using ForwardDiff
using LinearAlgebra
using Revise

export get_line_quadrature
export get_all_line_quadratures

function get_line_quadrature(ref_quadr, curve, s_domain, curveParams...; ref_domain=(-1,1), normalization=false)

    # Map the reference quadr points to s values
    ref_pts, ref_weights = ref_quadr
    x_lb, x_ub = ref_domain
    s_lb, s_ub = s_domain
    s_pts = @. (s_ub - s_lb)/(x_ub - x_lb)*(ref_pts - x_lb) + s_lb

    # Calculate the derivative and normal values at each of the s values
    dc_ds(s) = ForwardDiff.derivative(t->curve(t,curveParams...), s)
    tangents = @. dc_ds(s_pts)


    normals = [ [ tangents[i][2], -tangents[i][1] ] for i in 1:length(tangents) ]
    if normalization == true
        for i = 1:length(normals)
            normals[i] = normals[i] / norm(normals[i])
        end
    end

    # Calculate the evaluation points and adjusted weights
    line_pts = @. curve(s_pts, curveParams...)
    line_weights = @. ref_weights.*norm(dc_ds(s_pts))

    return line_pts, line_weights, normals
end

function get_all_line_quadratures(ref_quadr, curve, stop_points, curveParams; ref_domain=(-1,1))

end

end # module
