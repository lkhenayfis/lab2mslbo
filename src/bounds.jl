
"""
    fix_infinity_bounds!(v::Vector{Float64}, sub::Float64)

Change `Inf` upper bounds to `sub`, as `DualSDDP` doesn't allow unbounded variables
"""
function fix_infinity_bounds!(md::ModelData, sub_lower::Float64, sub_upper::Float64)
    for i in range(1, length(md.x_upper))
        md.x_upper[i] = isinf(md.x_upper[i]) ? sub_upper : md.x_upper[i]
    end

    for i in range(1, length(md.x_lower))
        md.x_lower[i] = isinf(md.x_lower[i]) ? sub_lower : md.x_lower[i]
    end
end

