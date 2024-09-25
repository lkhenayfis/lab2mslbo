using Random
using SDDPlab
using SDDP, JuMP

include("helpers.jl")

function build_mslbo(data_dir::String, seed::Integer = 1234)::DualSDDP.MSLBO

    Random.seed!(seed)
    aux, saa = build_sddp_model(data_dir)
    # loop em aux.nodes extraindo as matrizes
    # builda um MSLBO
end

"""
    build_sddp_model

Build a SDDP.PolicyGraph from which system matrices are extracted

# Argumentos

  - `data_dir::String`: diretorio de um caso do lab
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

function jump2matrices(m::JuMP.Model)::Vector{Matrix}

    x, full_A, c, lb, ub = split_elements(m)

    # cleans out slacks
    # this is necessary because dualsddp is built for equality constraints only, which would demand
    # creating slacks of slacks after building the JuMP Model
    # too much work, can be fixed in inputs and should be better addressed in the backend
    full_A, lb, ub = remove_variable_by_name(variables, "_SLACK", full_A, lb, ub)
    
    # outflow and net exchange are convenience variables meant to simplify outputs
    # the reason they are cut off from the problem is for being unbounded, which causes problems in
    # the dualsddp implementation
    # given they serve no actual modeling purpose, it's easier to just remove them
    full_A, lb, ub = remove_variable_by_name(variables, "_OUTFLOW", full_A, lb, ub)
    full_A, lb, ub = remove_variable_by_name(variables, "_NET_EXCHANGE", full_A, lb, ub)

    # cleans out dummy variables (added but not used in any constraint)
    variables, full_A, costs = remove_dummy_variables(variables, full_A, costs)

    bounds, affine = split_bounds_affine(full_A, lb, ub)

    # checks if any affine constraints

end

# AUXILIARES ---------------------------------------------------------------------------------------

function get_costs(m::JuMP.Model, variables::Vector{VariableRef})
    terms = m.ext[:sddp_node].stage_objective.terms
    costs = zeros(Float64, size(variables))
    for (v, coef) in terms
        index = findfirst(==(v), variables)
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

    return variables, costs, A, b_lower, b_upper
end

"""
    split_full_A(mat::Matrix)::Tuple{Matrix}

Split a technology matrix A into submatrix of bound constraints and general constraints
"""
function split_bounds_affine(technology::Matrix, lower::Vector, upper::Vector)::Tuple

    lines_of_one = [is_line_of_one(r) for r in eachrow(technology)]

    A_bounds = technology[lines_of_one, :]
    A_affine = technology[.!lines_of_one, :]

    l_bounds = lower[lines_of_one]
    l_affine = lower[.!lines_of_one]

    u_bounds = upper[lines_of_one]
    u_affine = upper[.!lines_of_one]

    return ((A_bounds, l_bounds, u_bounds), (A_affine, l_affine, u_affine))
end

function is_line_of_one(row)::Bool
    num_not_zero = sum(row .!= 0.0)
    equals_one = sum(row .== 1.0)

    return (num_not_zero == 1) & (equals_one == 1)
end

