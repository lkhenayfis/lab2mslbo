
"""
    fix_unbounded(bound_mat::Matrix, lb::Vector{Float64}, ub::Vector{Float64})

Add upper bound at Infinity for any unbounded variables in the problem

# Arguments

All arguments are, in order, the first returned `Tuple` from `split_bounds_affine`

  - `bound_mat::Matrix` section of technology matrix regarding box constraints
  - `lb::Vector{Float64}` section of lower bounds regarding box constraints
  - `ub::Vector{Float64}` section of upper bounds regarding box constraints

# Return

A `Tuple` of it's arguments, in the same order, after inclusion of new bounds
"""
function fix_unbounded(bound_mat::Matrix, lb::Vector{Float64}, ub::Vector{Float64})::Tuple
    Nvars = size(bound_mat)[2]
    unbounded = mapslices(x -> sum(x), bound_mat, dims=[1])[1,:]
    unbounded = findall(unbounded .== 0.0)

    for u in unbounded
        newrow = [i == u ? 1 : 0 for i in 1:Nvars]
        bound_mat = cat(bound_mat, newrow'; dims=1)
    end

    lb = vcat(lb, repeat([0], length(unbounded)))
    ub = vcat(ub, repeat([Inf], length(unbounded)))

    return bound_mat, lb, ub
end

"""
    fix_infinity_bounds!(v::Vector{Float64}, sub::Float64)

Change `Inf` upper bounds to `sub`, as `DualSDDP` doesn't allow unbounded variables
"""
function fix_infinity_bounds!(v::Vector{Float64}, sub::Float64)
    for i in range(1, length(v))
        v[i] = isinf(v[i]) ? sub : v[i]
    end
end

"""
    sort_by_selector_matrix(sel_mat::Matrix, v::Vector{Float64})::Vector{Float64}

Return a reordering of vector `v` based on a selector matrix `sel_mat`

Following the discussion in `split_bounds_affine`, there is still one additional complication
regarding box constraints. It's section of the technology matrix isn't ordered to be diagonal, thus
the vectors of lower and upper bounds are misaligned with the vector of variables. Given in
`DualSDDP` we never declare variables directly, only system matrices from which variables are
inferred, consistency in order of variables must be mantained between affine and box. This function
exists to reorder lower and upper bounds accordingly.

# Arguments

  - `sel_mat::Matrix` a selection matrix, i.e. the section of box contraints from the technology
    matrix
  - `v::Vector{Float64}` a vector to be reordered
"""
function sort_by_selector_matrix(sel_mat::Matrix, v::Vector{Float64})::Vector{Float64}
    order_sel = Vector{Int}()

    for row in eachrow(sel_mat)
        ind = findfirst(row .== 1.0)
        push!(order_sel, ind)
    end

    order = sortperm(order_sel)
    sorted_v = v[order]

    return sorted_v
end

"""
    split_bounds(vars, bound_mat, lb, ub)

Split lower and upper bounds vectors into ones for state and control variables separately

# Arguments

  - `vars::Vector{VariableRef}` subproblem variables
  - `bound_mat::Matrix` technology matrix section of box constraints as returned in
    `split_bounds_affine`
  - `lb::Vector{Float64}` section of lower bounds regarding box constraints
  - `ub::Vector{Float64}` section of upper bounds regarding box constraints

# Returns

A `Tuple` of two `Tuples` each containing lower and upper bounds (in this order) for state and
control variables (in this order as well)
"""
function split_bounds(vars::Vector{VariableRef}, bound_mat::Matrix,
    lb::Vector{Float64}, ub::Vector{Float64})::Tuple

    state_t, state_t_1, control = get_state_control_indexes(vars)

    sorted_lb = sort_by_selector_matrix(bound_mat, lb)
    sorted_ub = sort_by_selector_matrix(bound_mat, ub)

    lb_states = vec(sorted_lb[state_t,:])
    ub_states = vec(sorted_ub[state_t,:])

    lb_control = vec(sorted_lb[control,:])
    ub_control = vec(sorted_ub[control,:])

    return (lb_states, ub_states), (lb_control, ub_control)
end
