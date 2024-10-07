
mutable struct ModelData
    x::Vector{VariableRef}
    c::Vector{Float64}
    A::Matrix{Float64}
    b_lower::Vector{Float64}
    b_upper::Vector{Float64}
    x_lower::Vector{Float64}
    x_upper::Vector{Float64}
end

function ModelData(m::JuMP.Model)
    split = split_elements(m)
    ModelData(split...)
end

# HELPERS ------------------------------------------------------------------------------------------

"""
    get_costs(m::JuMP.Model, vars::Vector{VariableRef})::Vector

Extract vector of costs from a `JuMP.Model` built for a `SDDP.PolicyGraph`
"""
function get_costs(m::JuMP.Model, vars::Vector{VariableRef})::Vector
    terms = m.ext[:sddp_node].stage_objective.terms
    costs = zeros(Float64, size(vars))
    for (v, coef) in terms
        index = findfirst(==(v), vars)
        costs[index] = coef
    end

    return costs
end

"""
    split_elements(m::JuMP.Model)::Tuple

Split a `JuMP.Model` into it's composing elements: costs, constraints matrix and bounds

# Arguments

  - `m::JuMP.Model` model built for a `SDDP.PolicyGraph`, as returned in nodes of `build_sddp_model`

# Returns

A `Tuple` containing, in order: vector of subproblem variables (`Vector{VariableRef}`), technology
matrix, vector of costs, vector of lower and upper bounds on the constraints, as in
described in `JuMP.lp_matrix_data`
"""
function split_elements(m::JuMP.Model)::Tuple
    
    md = lp_matrix_data(m)

    # always removes epigraphical variable added last by SDDP.jl
    variables = md.variables[1:end-1]
    costs = get_costs(m, variables)
    A = Matrix(md.A)[:, 1:end-1]
    b_lower = md.b_lower
    b_upper = md.b_upper
    x_lower = md.x_lower[1:end-1]
    x_upper = md.x_upper[1:end-1]

    return variables, costs, A, b_lower, b_upper, x_lower, x_upper
end

# METHODS ------------------------------------------------------------------------------------------

"""
    get_state_control_indexes(md::ModelData)::Tuple

Return `Tuple` of three vectors: position of current state, last state and control variables in the
vector of variables extracted from the `JuMP.Model`
"""
function get_state_control_indexes(md::ModelData)::Tuple
    state_t = get_variable_index(md.x, "_out")
    state_t_1 = get_variable_index(md.x, "_in")

    control = Vector{Int}()
    for i in range(1, length(md.x))
        !((i in state_t) | (i in state_t_1)) ? push!(control, i) : nothing
    end

    return state_t, state_t_1, control
end

"""
    split_constraint_matrices(md::ModelData)::Tuple

Return `Tuple` of matrices `A`, `B` and `T` as specified in `DualSDDP`

# Arguments

  - `vars::Vector{VariableRef}` subproblem variables
  - `tech_mat::Matrix` technology matrix OF THE AFFINE SECTION, as returned from
    `split_bounds_affine`
"""
function split_constraint_matrices(md::ModelData)
    state_t, state_t_1, control = get_state_control_indexes(md)

    A = md.A[:, state_t]
    B = md.A[:, state_t_1]
    T = md.A[:, control]

    return A, B, T
end

function get_equality_inequality_indexes(md::ModelData)::Tuple
    equality = findall(md.b_lower .== md.b_upper)
    inequality = findall(md.b_lower .!= md.b_upper)

    return equality, inequality
end
