
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

# OUTPUT GENERATION --------------------------------------------------------------------------------

function export_primal_with_ub_convergence(
    niters::Int64,
    primal_lbs::Vector{Float64},
    primal_times::Vector{Float64},
    rec_ubs::Vector{Tuple{Int64,Float64}},
    rec_times::Vector{Float64},
    writer::Function,
    extension::String;
    output_path_without_extension::String = "./convergence",
)
    dense_ubs = fill(NaN, niters)
    for ub_pair in rec_ubs
        dense_ubs[ub_pair[1]] = ub_pair[2]
    end
    dense_rec_times = fill(NaN, niters)
    for (ub_pair, time) in zip(rec_ubs, rec_times)
        dense_rec_times[ub_pair[1]] = time
    end

    df = DataFrame(;
        iteration = 1:niters,
        lower_bound = primal_lbs,
        simulation = fill(NaN, niters),
        upper_bound = dense_ubs,
        primal_time = primal_times,
        upper_bound_time = dense_rec_times,
        time = primal_times + replace(dense_rec_times, NaN => 0.0),
    )

    return writer(output_path_without_extension * extension, df)
end

function export_primal_with_dual_ub_convergence(
    niters::Int64,
    primal_lbs::Vector{Float64},
    primal_times::Vector{Float64},
    dual_ubs::Vector{Float64},
    dual_times::Vector{Float64},
    writer::Function,
    extension::String;
    output_path_without_extension::String = "./convergence",
)
    df = DataFrame(;
        iteration = 1:niters,
        lower_bound = primal_lbs,
        simulation = fill(NaN, niters),
        upper_bound = dual_ubs,
        primal_time = primal_times,
        upper_bound_time = dual_times,
        time = primal_times + dual_times,
    )

    return writer(output_path_without_extension * extension, df)
end

function export_problem_child_convergence(
    niters::Int64,
    io_lbs::Vector{Float64},
    io_ubs::Vector{Float64},
    io_times::Vector{Float64},
    writer::Function,
    extension::String;
    output_path_without_extension::String = "./convergence",
)
    df = DataFrame(;
        iteration = 1:niters,
        lower_bound = io_lbs,
        simulation = fill(NaN, niters),
        upper_bound = io_ubs,
        time = io_times,
    )

    return writer(output_path_without_extension * extension, df)
end

function translate_cuts_states(filename::String, vertex_name_parser::Function)
    policy = JSON.parsefile(filename; use_mmap = false)
    keys_to_translate = ["state", "coefficients"]
    for node_info in policy
        for cut_info in node_info["single_cuts"]
            for key in keys_to_translate
                new_info = Dict{String,Float64}([])
                for (name, value) in cut_info[key]
                    new_info[vertex_name_parser(name)] = value
                end
                cut_info[key] = new_info
            end
        end
    end
    open(filename, "w") do f
        JSON.print(f, policy)
    end
    return nothing
end

function primal_cuts_json_to_table(
    pb::DualSDDP.MSSP,
    writer::Function,
    extension::String;
    output_path_without_extension::String = "./cuts",
) end

function dual_cuts_json_to_table(
    pb::DualSDDP.MSSP,
    writer::Function,
    extension::String;
    output_path_without_extension::String = "./cuts",
) end