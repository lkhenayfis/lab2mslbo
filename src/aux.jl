
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