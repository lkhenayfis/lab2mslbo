
"""
    build_mslbo(data_dir::String; seed::Integer = 1234)::DualSDDP.MSLBO

Build a `DualSDDP.MSLBO` object from a `SDDPlab` input data directory
"""
function build_mslbo(data_dir::String; seed::Integer = 1234)::DualSDDP.MSLBO

    Random.seed!(seed)
    aux, saa = build_sddp_model(data_dir)
    
    lbos = Vector{DualSDDP.SimpleLBO}()
    for i in 1:length(aux.nodes)
        A, B, T, c, ds, sb, cb = jump2matrices(aux.nodes[i].subproblem, saa[i])
        probs = repeat([1], length(ds)) ./ length(ds)
        lbo = DualSDDP.SimpleLBO(A, B, T, c, ds, sb[2], cb[2], 0, 1e5, 1e5, probs)
        push!(lbos, lbo)
    end

    return DualSDDP.build(lbos)
end

"""
    build_sddp_model

Build a SDDP.PolicyGraph from which system matrices are extracted

# Arguments

  - `data_dir::String`: path to a `SDDPlab` case data directory

# Returns

A `Tuple` containing two elements: the first is an `SDDP.PolicyGraph` containing the full model;
second is the Sample Average Aproximation built, i.e., noises for every stage, branch and random
element (in that order) in a `Vector{Vector{Vector}}`
"""
function build_sddp_model(data_dir::String)::Tuple
    original_wd = pwd()

    cd(data_dir)

    entrypoint = SDDPlab.Inputs.Entrypoint("main.jsonc", CompositeException())
    path = SDDPlab.Inputs.get_path(entrypoint)
    files = SDDPlab.Inputs.get_files(entrypoint)

    cd(original_wd)

    model = SDDPlab.Tasks.__build_model(files)

    scenarios = SDDPlab.Inputs.get_scenarios(files)
    num_stages = SDDPlab.Inputs.get_number_of_stages(SDDPlab.Inputs.get_algorithm(files))

    SAA = SDDPlab.Scenarios.generate_saa(scenarios, num_stages)

    return model, SAA
end

"""
    jump2matrices(m::JuMP.Model, node_noises::Vector{Vector{Float64}})

Convert a JuMP model and noise sample space to MSLBO compliant matrices

# Arguments

  - `m::JuMP.Model` the model of a given node in a `SDDP.PolicyGraph`, as built by 
      `build_sddp_model`
  - `node_noises::Vector{Vector{Float64}}` sample space of noises in a given node; this is an 
      element of the SAA returned by `build_sddp_model`

# Returns

A `Tuple` containing, in order: matrices `A`, `B` and `T` as declared in `DualSDDP`; vector of
costs, vector of vector of `d`s as declared in `DualSDDP`, tuple of states' lower and upper
bounds and tuple of controls' lower and upper bounds, all vectors 
"""
function jump2matrices(m::JuMP.Model, node_noises::Vector{Vector{Float64}})::Tuple

    x, full_A, c, lb, ub = split_elements(m)
    x, full_A, c, lb, ub = sanitize_variables(x, full_A, c, lb, ub)

    x, full_A, c, ds, sp_c = fix_ω(x, full_A, c, ub, node_noises)

    bounds, affine = split_bounds_affine(full_A, lb, ub, ds, sp_c)
    bounds = fix_unbounded(bounds...)
    
    A, B, T = split_constraint_matrices(x, affine[1])

    states_bounds, control_bounds = split_bounds(x, bounds...)

    return A, B, T, c, ds, states_bounds, control_bounds
end

# AUXILIARES ---------------------------------------------------------------------------------------

function get_costs(m::JuMP.Model, vars::Vector{VariableRef})
    terms = m.ext[:sddp_node].stage_objective.terms
    costs = zeros(Float64, size(vars))
    for (v, coef) in terms
        index = findfirst(==(v), vars)
        costs[index] = coef
    end

    return costs
end

function split_elements(m::JuMP.Model)
    
    md = lp_matrix_data(m)

    # always removes epigraphical variable added last by SDDP.jl
    variables = md.variables[1:end-1]
    costs = get_costs(m, variables)
    A = Matrix(md.A)[:, 1:end-1]
    b_lower = md.b_lower
    b_upper = md.b_upper

    return variables, A, costs, b_lower, b_upper
end

function sanitize_variables(vars, tech_mat, costs, b_lower, b_upper)
    # cleans out slacks
    # this is necessary because dualsddp is built for equality constraints only, which would demand
    # creating slacks of slacks after building the JuMP Model
    # too much work, can be fixed in inputs and should be better addressed in the backend
    tech_mat, b_lower, b_upper = remove_variable_by_name(vars, "_SLACK", tech_mat, b_lower, b_upper)
    
    # outflow and net exchange are convenience variables meant to simplify outputs
    # the reason they are cut off from the problem is for being unbounded, which causes problems in
    # the dualsddp implementation
    # given they serve no actual modeling purpose, it's easier to just remove them
    #tech_mat, b_lower, b_upper = remove_variable_by_name(vars, "OUTFLOW", tech_mat, b_lower, b_upper) # if dropped breaks hydro balance
    tech_mat, b_lower, b_upper = remove_variable_by_name(vars, "NET_EXCHANGE", tech_mat, b_lower, b_upper)

    # cleans out dummy variables (added but not used in any constraint)
    vars, tech_mat, costs = remove_dummy_variables(vars, tech_mat, costs)

    return vars, tech_mat, costs, b_lower, b_upper
end

function is_line_of_one(row)::Bool
    num_not_zero = sum(row .!= 0.0)
    equals_one = sum(row .== 1.0)

    return (num_not_zero == 1) & (equals_one == 1)
end

"""
    split_full_A(mat::Matrix)::Tuple{Matrix}

Split a technology matrix A into submatrix of bound constraints and general constraints
"""
function split_bounds_affine(tech_mat::Matrix, b_lower::Vector, b_upper::Vector,
    ds::Vector{Vector{Float64}}, except::Vector{Bool})::Tuple

    lines_of_one = [is_line_of_one(r) for r in eachrow(tech_mat)]
    lines_of_one[except] .= false

    A_bounds = tech_mat[lines_of_one, :]
    A_affine = tech_mat[.!lines_of_one, :]

    l_bounds = b_lower[lines_of_one]
    l_affine = b_lower[.!lines_of_one]

    u_bounds = b_upper[lines_of_one]
    ds = [d[.!lines_of_one] for d in ds]

    return ((A_bounds, l_bounds, u_bounds), (A_affine, l_affine, ds))
end

function fix_unbounded(bound_mat::Matrix, lb::Vector{Float64}, ub::Vector{Float64})
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

function get_state_control_indexes(vars::Vector{VariableRef})
    state_t = get_variable_index(vars, "_in")
    state_t_1 = get_variable_index(vars, "_out")

    control = Vector{Int}()
    for i in range(1, length(vars))
        !((i in state_t) | (i in state_t_1)) ? push!(control, i) : nothing
    end

    return state_t, state_t_1, control
end

function split_constraint_matrices(vars::Vector{VariableRef}, tech_mat::Matrix)
    state_t, state_t_1, control = get_state_control_indexes(vars)

    A = tech_mat[:, state_t]
    B = tech_mat[:, state_t_1]
    T = tech_mat[:, control]

    return A, B, T
end

function sort_by_selector_matrix(sel_mat::Matrix, v::Vector{Float64})
    order_sel = Vector{Int}()

    for row in eachrow(sel_mat)
        ind = findfirst(row .== 1.0)
        push!(order_sel, ind)
    end

    order = sortperm(order_sel)
    sorted_v = v[order]

    return sorted_v
end

function split_bounds(vars::Vector{VariableRef}, bound_mat::Matrix,
    lb::Vector{Float64}, ub::Vector{Float64})

    state_t, state_t_1, control = get_state_control_indexes(vars)
    states = vcat(state_t, state_t_1)

    sorted_lb = sort_by_selector_matrix(bound_mat, lb)
    sorted_ub = sort_by_selector_matrix(bound_mat, ub)

    lb_states = vec(sorted_lb[states,:])
    ub_states = vec(sorted_ub[states,:])

    lb_control = vec(sorted_lb[control,:])
    ub_control = vec(sorted_ub[control,:])

    return (lb_states, ub_states), (lb_control, ub_control)
end
