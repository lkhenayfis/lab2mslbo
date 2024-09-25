# MISC ---------------------------------------------------------------------------------------------

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

# SANITAZING SLACK VARIABLES -----------------------------------------------------------------------

"""
    parse_slack_variables

Find slack variables and remove them from the system matrices
"""
function parse_slack_variables(variables::Vector{VariableRef},
    A::Matrix, lower::Vector, upper::Vector)

    pattern = "_SLACK"
    indexes = get_variable_index(variables, pattern)
    slack_variables = variables[indexes]

    if length(indexes) >= 1
        warn_slacks(slack_variables)
        bool_drop = .!find_non_zero(indexes, A)
        A = A[bool_drop, :]
        lower = lower[bool_drop]
        upper = upper[bool_drop]
    end

    return A, lower, upper
end

"""
    find_non_zero

Return vector of bools pointing in which lines of `A` at least one column in `i` has non zero value
"""
function find_non_zero(i::Vector{Int}, A::Matrix)
    has_nonzero = Vector{Bool}()
    for row in eachrow(A)
        push!(has_nonzero, sum(row[i]) > 0)
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
    parse_dummy_variables

Find dummy variables (not used in any constraint) and remove them from the system matrices
"""
function parse_dummy_variables(variables::Vector{VariableRef}, A::Matrix, costs::Vector)

    indexes = find_dummies(A)
    dummy_variables = variables[indexes]

    if length(indexes) >= 1
        warn_dummies(dummy_variables)
        variables = variables[.!indexes]
        A = A[:, .!indexes]
        costs = costs[.!indexes]
    end

    return variables, A, costs
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

# LOWER/UPPER --------------------------------------------------------------------------------------

function check_equal_affine_bounds(lower::Vector, upper::Vector)
    equals = [l == b for (l, b) in zip(lower, upper)]
    if sum(equals) != length(equals)
        @warn "Some affine constraints are not of equality type -- revise inputs"
    end
end

function check_zero_lower_bound(lower::Vector, variables::Vector{VariableRef})
    is_zero = [l == 0 for l in lower]
    if sum(is_zero) != length(is_zero)
        warn_msg = "Variables " * paste_vector(string.(variables[.!is_zero]))
        warn_msg *= " have non-zero lower bound -- revise inputs"
        @warn warn_msg
    end
end