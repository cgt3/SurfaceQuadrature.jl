module SurfaceQuadrature

using ForwardDiff
using Revise

export get_line_quadrature


function get_line_quadrature(ref_quadr, curve, s_domain, curveParams...; ref_domain=(-1,1))
    x_lb, x_ub = ref_domain
    s_lb, s_ub = s_domain

    ref_pts, ref_weights = ref_quadr
    # TODO: should return the mapped s values or the (x,y...) points associated with the curve?
    line_pts = @. ref_pts + (s_ub - s_lb)/(x_ub - x_lb)*(ref_pts - x_lb) + s_lb

    dc_ds(s) = ForwardDiff.derivative(t->curve(t,curveParams...), s)
    line_weights = ref_weights.*dc_ds(line_points)

    return (line_pts, line_weights)
end

end
