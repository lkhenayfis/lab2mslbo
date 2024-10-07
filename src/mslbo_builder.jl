
"""
    build_mslbo(data_dir::String; seed::Integer = 1234)::DualSDDP.MSLBO

Build a `DualSDDP.MSLBO` object from a `SDDPlab` input data directory

# Arguments

  - `data_dir::String` directory containig `SDDPlab` compliant input files
  - `problem_ub::Vector{Float64}`: Upper bound on the value of the problem at each stage
  - `α::Vector{Float64}`: Upper bound on the Lipschitz constant at each stage
"""
function build_mslbo(data_dir::String;
    fun_problem_ub::Function = default_guess,
    fun_α::Function = default_guess,
    seed::Integer = 1234)::DualSDDP.MSLBO

    files = read_lab_inputs(data_dir)
    problem_ub = fun_problem_ub(files)
    α = fun_α(files)

    Random.seed!(seed)
    aux, saa = build_sddp_model(files)
    
    lbos = Vector{DualSDDP.SimpleLBO}()
    for i in 1:length(aux.nodes)
        A, B, T, c, ds, sb, cb = jump2matrices(aux.nodes[i].subproblem, saa[i])
        probs = repeat([1], length(ds)) ./ length(ds)
        lbo = DualSDDP.SimpleLBO(A, B, T, c, ds, sb[2], cb[2], 0, problem_ub, α, probs)
        push!(lbos, lbo)
    end

    return DualSDDP.build(lbos)
end

function default_guess(files::Vector{SDDPlab.Inputs.InputModule})
    algo = SDDPlab.Inputs.get_algorithm(files)
    num_stages = SDDPlab.Inputs.get_number_of_stages(algo)

    buses = SDDPlab.System.get_buses_entities(SDDPlab.Inputs.get_system(files))
    scens = SDDPlab.Inputs.get_scenarios(files)

    guess = 0.0
    for stage in range(1, num_stages)
        for bus in buses
            guess += bus.deficit_cost * SDDPlab.Scenarios.get_load(bus.id, stage, scens)
        end
    end

    return guess
end

function read_lab_inputs(data_dir::String)::Vector{SDDPlab.Inputs.InputModule}
    original_wd = pwd()

    cd(data_dir)

    entrypoint = SDDPlab.Inputs.Entrypoint("main.jsonc", CompositeException())
    path = SDDPlab.Inputs.get_path(entrypoint)
    files = SDDPlab.Inputs.get_files(entrypoint)

    cd(original_wd)

    return files
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
function build_sddp_model(files::Vector{SDDPlab.Inputs.InputModule})::Tuple
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
    c = reduce_cost(x, c)

    states_bounds, control_bounds = split_bounds(x, bounds...)
    fix_infinity_bounds!(control_bounds[2], 1e4)

    return A, B, T, c, affine[3], states_bounds, control_bounds
end

# AUXILIARES ---------------------------------------------------------------------------------------

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

    return variables, A, costs, b_lower, b_upper
end

"""
    sanitize_variables(vars, tech_mat, costs, b_lower, b_upper)

Clean some variables to fit within the `DualSDDP` problem specification

Specifically, this function removes variables and constraints in which they are included, otherwise
errors would be induced.

It first removes slack variables, which don't fit within `DualSDDP`'s problem specification directly
(demands some manipulation and inclusion of slacks of slacks, after building the JuMP.Models which
is not worth the trouble). Then the `NET_EXCHANGE` variable is removed, as it's simply a convenience
variable and should be converted into an AffineExpression in JuMP eventually. Finally are removed
what are called 'dummy variables', i.e. variables which don't appear in any constraint. This can
happen after earlier removals and also is the case of `LOAD`, which is included as a decision
variable in preparation for it being stochastic.

# Arguments

All arguments of this function are the elements returned by `split_elements`, in the same order:

  - `vars::Vector{VariableRef}` problem variables
  - `tech_mat::Matrix` technology matrix
  - `costs::Vector{Float64}` vector of costs
  - `b_lower::Vector{Float64}` constraints lower bounds
  - `b_upper::Vector{Float64}` constraints upper bounds

# Returns

A `Tuple` of it's arguments, in the same order, after removal of variables
"""
function sanitize_variables(vars::Vector{VariableRef}, tech_mat::Matrix,
    costs::Vector{Float64}, b_lower::Vector{Float64}, b_upper::Vector{Float64})::Tuple
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

"""
    is_line_of_one(row)::Bool

Check if a matrix line `row` has only one element and it is equal to 1
"""
function is_line_of_one(row)::Bool
    num_not_zero = sum(row .!= 0.0)
    equals_one = sum(row .== 1.0)

    return (num_not_zero == 1) & (equals_one == 1)
end

"""
    split_bounds_affine

Split composing elements into the `affine` and `bounds` section

Usually `JuMP.lp_matrix_data` returns box constraints in a separate element from the upper and lower
bounds on constraints. However, `SDDPlab` declares variables in a way which doesn't allow for this
representation, thus mixing into the model what are actual affine constraints and simple box ones.
This function splits the tecnology matrix, lower and upper bounds on constraints into these two
sections

The final two arguments benefit from extra details. `ds` is a Vector{Vector} where inner vectors are
but the rhs of equality constraints, but with the term corresponding to the stochastic process
constraint modified to a noise value. This collection of vectors represents all possible variations
of rhs in a given stage. `except` is a convenience argument, needed because in the simplest cases
the stochastic process constraint appears to be a bounds one, as it only concerns one variable with
corresponding coefficient 1.

# Arguments

  - `tech_mat::Matrix` technology matrix
  - `b_lower::Vector{Float64}` constraints lower bounds
  - `b_upper::Vector{Float64}` constraints upper bounds
  - `ds::Vector{Vector{Float64}}` variations of the rhs in equality constraints
  - `except::Vector{Bool}` a vector indicating which lines of `tech_mat` should be excluded from
    inference of wether it is a bounds or affine constraint (used directly as affine)

# Return

A `Tuple` containing two `Tuple`s: the first contains the box constraints section of the technology
matrix, lower and upper; the second is the safe, but for affine constraints, except for the final
element: instead of upper, it is the collection `ds` subseted in the positions of affine constraints
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

"""
    split_constraint_matrices(vars::Vector{VariableRef}, tech_mat::Matrix)::Tuple

Return `Tuple` of matrices `A`, `B` and `T` as specified in `DualSDDP`

# Arguments

  - `vars::Vector{VariableRef}` subproblem variables
  - `tech_mat::Matrix` technology matrix OF THE AFFINE SECTION, as returned from
    `split_bounds_affine`
"""
function split_constraint_matrices(vars::Vector{VariableRef}, tech_mat::Matrix)
    state_t, state_t_1, control = get_state_control_indexes(vars)

    A = tech_mat[:, state_t]
    B = tech_mat[:, state_t_1]
    T = tech_mat[:, control]

    return A, B, T
end

function reduce_cost(vars::Vector{VariableRef}, costs::Vector{Float64})
    state_t, state_t_1, control = get_state_control_indexes(vars)
    costs = costs[control]
    return costs
end
