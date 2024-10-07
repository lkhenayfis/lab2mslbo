
"""
    get_variable_index(v::Vector{VariableRef}, p::String)

Look for pattern `p` in vector `v`, returning indexes of matched elements
"""
function get_variable_index(v::Vector{VariableRef}, p::String)
    v_string = string.(v)
    matches = match.(Regex(p), v_string)

    index_matches = Vector{Int}()
    for (i, m) in enumerate(matches)
        !isnothing(m) ? push!(index_matches, i) : nothing
    end

    return index_matches
end

"""
    get_state_control_indexes(vars::Vector{VariableRef})::Tuple

Return `Tuple` of three vectors: position of current state, last state and control variables in the
vector of variables extracted from the `JuMP.Model`
"""
function get_state_control_indexes(vars::Vector{VariableRef})::Tuple
    state_t = get_variable_index(vars, "_out")
    state_t_1 = get_variable_index(vars, "_in")

    control = Vector{Int}()
    for i in range(1, length(vars))
        !((i in state_t) | (i in state_t_1)) ? push!(control, i) : nothing
    end

    return state_t, state_t_1, control
end

"""
    paste_vector

Take a vector of strings and return a single string `"(v_1, v_2, ..., v_n)"`
"""
function paste_vector(v::Vector{String})
    pasted = "(" * v[1]
    enum = Iterators.drop(enumerate(v), 1)
    
    for (ind, elem) in enum
        pasted *= (", " * v[ind])
    end
    pasted *= ")"

    return pasted
end

# SANITAZING VARIABLE BY NAME ----------------------------------------------------------------------

"""
    remove_variable_by_name

Find slack variables and remove them from the system matrices
"""
function remove_variable_by_name(vars::Vector{VariableRef}, pattern::String,
    tech_mat::Matrix, b_lower::Vector, b_upper::Vector)

    indexes = get_variable_index(vars, pattern)
    slack_variables = vars[indexes]

    if length(indexes) >= 1
        warn_slacks(slack_variables)
        bool_drop = .!find_non_zero(indexes, tech_mat)
        tech_mat = tech_mat[bool_drop, :]
        b_lower = b_lower[bool_drop]
        b_upper = b_upper[bool_drop]
    end

    return tech_mat, b_lower, b_upper
end

"""
    find_non_zero

Return vector of bools pointing in which lines of `A` at least one column in `i` has non zero value
"""
function find_non_zero(i::Vector{Int}, A::Matrix)
    has_nonzero = Vector{Bool}()
    for row in eachrow(A)
        push!(has_nonzero, sum(row[i] .!= 0) > 0)
    end

    return has_nonzero
end

function warn_slacks(sv::Union{VariableRef,Vector{VariableRef}})
    warn_msg = "Found the following slack variables: " * paste_vector(string.(sv))
    warn_msg *= "\n these and their associated constraints will be removed -- revise the inputs"
    @warn warn_msg
end

# SANITIZING DUMMY VARIABLES -----------------------------------------------------------------------

"""
    remove_dummy_variables

Find dummy variables (not used in any constraint) and remove them from the system matrices
"""
function remove_dummy_variables!(md::ModelData)

    indexes = find_dummies(md.A)
    dummy_variables = md.x[indexes]

    if length(indexes) >= 1
        warn_dummies(dummy_variables)
        md.x = md.x[.!indexes]
        md.A = md.A[:, .!indexes]
        md.c = md.c[.!indexes]
        md.x_lower = md.x_lower[.!indexes]
        md.x_upper = md.x_upper[.!indexes]
    end
end

"""
    find_dummies

Return vector of bools indicating which columns of `A` are fully zero
"""
function find_dummies(A::Matrix)
    is_dummy = Vector{Bool}()
    for col in eachcol(A)
        push!(is_dummy, all([elem == 0 for elem in col]))
    end

    return is_dummy
end

function warn_dummies(dv::Union{VariableRef,Vector{VariableRef}})
    warn_msg = "Found the following dummy variables: " * paste_vector(string.(dv))
    warn_msg *= "\n these and their associated constraints will be removed -- revise the inputs"
    @warn warn_msg
end

# NOISES -------------------------------------------------------------------------------------------

function fix_ω(vars, tech_mat, costs, b_upper, node_noises)
    # zeroes out any coefficients in tech_mat related to noise terms
    ω_indices = get_variable_index(vars, "ω_")
    sp_constraints = find_non_zero(ω_indices, tech_mat)
    tech_mat[:, ω_indices] .= 0

    # rodar um remove_dummy_variables para limpar os ω_X
    vars, tech_mat, costs = remove_dummy_variables(vars, tech_mat, costs)

    B = length(node_noises)
    ds = [copy(b_upper) for i in range(1, B)]
    for (i, d) in enumerate(ds)
        d[sp_constraints] .= node_noises[i]
    end

    return vars, tech_mat, costs, ds, sp_constraints
end

# LOWER/UPPER --------------------------------------------------------------------------------------

function check_equal_affine_bounds(lower::Vector, upper::Vector)
    equals = [l == b for (l, b) in zip(lower, upper)]
    if sum(equals) != length(equals)
        @warn "Some affine constraints are not of equality type -- revise inputs"
    end
end

function check_zero_lower_bound(lower::Vector, vars::Vector{VariableRef})
    is_zero = [l == 0 for l in lower]
    if sum(is_zero) != length(is_zero)
        warn_msg = "Variables " * paste_vector(string.(vars[.!is_zero]))
        warn_msg *= " have non-zero lower bound -- revise inputs"
        @warn warn_msg
    end
end