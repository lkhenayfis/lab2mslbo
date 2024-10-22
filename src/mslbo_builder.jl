
"""
    build_mslbo(data_dir::String; seed::Integer = 1234)::DualSDDP.MSLBO

Build a `DualSDDP.MSLBO` object from a `SDDPlab` input data directory

# Arguments

  - `data_dir::String` directory containig `SDDPlab` compliant input files
  - `problem_ub::Vector{Float64}`: Upper bound on the value of the problem at each stage
  - `α::Vector{Float64}`: Upper bound on the Lipschitz constant at each stage
"""
function build_mslbo(data_dir::String;
    fun_problem_ub::Function=default_guess,
    fun_α::Function=default_guess)::Tuple{DualSDDP.MSLBO,LabData,Any,Function,String}

    files = read_lab_inputs(data_dir)
    seed = get_seed(files)
    problem_ub = fun_problem_ub(files)
    α = fun_α(files)

    Random.seed!(seed)
    aux, saa = build_sddp_model(files)
    initial_states = get_initial_states(files)

    lbos = Vector{DualSDDP.SimpleLBO}()
    for i in 1:length(aux.nodes)
        slbod = SimpleLBOData(aux.nodes[i].subproblem, saa[i], 0.0, problem_ub, α)
        lbo = DualSDDP.SimpleLBO(
            slbod.A_eq,
            slbod.B_eq,
            slbod.T_eq,
            slbod.c,
            slbod.d_eq,
            slbod.Ux,
            slbod.Uy,
            slbod.lb,
            slbod.ub,
            slbod.α,
            slbod.prob
        )
        push!(lbos, lbo)
    end

    num_stages = get_num_stages(files)
    risk_parameters = get_risk_measure_parameters(files)
    num_iterations = get_num_iterations(files)
    output_path = get_output_path(files)
    data = LabData(initial_states, seed, num_stages, num_iterations, output_path, risk_parameters[1], risk_parameters[2])

    return DualSDDP.build(lbos), data, get_solver(files), get_writer(files), get_extension(files)
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

    e = CompositeException()
    entrypoint = SDDPlab.Inputs.Entrypoint("main.jsonc", e)

    if length(e) > 0
        for exc in e
            println(exc)
        end
    end

    path = SDDPlab.Inputs.get_path(entrypoint)
    files = SDDPlab.Inputs.get_files(entrypoint)

    cd(original_wd)

    return files
end

function get_initial_states(files::Vector{SDDPlab.Inputs.InputModule})
    hydros = SDDPlab.System.get_hydros_entities(SDDPlab.Inputs.get_system(files))
    initial_states = [h.initial_storage for h in hydros]

    return initial_states
end


function get_seed(files::Vector{SDDPlab.Inputs.InputModule})
    return SDDPlab.Inputs.get_scenarios(files).seed
end

function get_num_stages(files::Vector{SDDPlab.Inputs.InputModule})
    return SDDPlab.Inputs.get_number_of_stages(SDDPlab.Inputs.get_algorithm(files))
end

function get_solver(files::Vector{SDDPlab.Inputs.InputModule})
    resources = SDDPlab.Inputs.get_resources(files)
    return SDDPlab.Inputs.generate_optimizer(resources.solver)
end

function get_num_iterations(files::Vector{SDDPlab.Inputs.InputModule})
    return SDDPlab.get_tasks(files)[2].convergence.max_iterations
end

function get_writer(files::Vector{SDDPlab.Inputs.InputModule})::Function
    result_format = SDDPlab.get_tasks(files)[2].results.format
    return SDDPlab.get_writer(result_format)
end

function get_extension(files::Vector{SDDPlab.Inputs.InputModule})::String
    result_format = SDDPlab.get_tasks(files)[2].results.format
    return SDDPlab.get_extension(result_format)
end

function get_output_path(files::Vector{SDDPlab.Inputs.InputModule})::String
    return SDDPlab.get_tasks(files)[2].results.path
end


function get_risk_measure_parameters(files::Vector{SDDPlab.Inputs.InputModule})
    risk_obj = SDDPlab.get_tasks(files)[2].risk_measure
    if nameof(typeof(risk_obj)) == :Expectation
        params = [1.0, 1.0]
    elseif nameof(typeof(risk_obj)) == :AVaR
        params = [risk_obj.alpha, 1.0]
    elseif nameof(typeof(risk_obj)) == :CVaR
        params = [risk_obj.alpha, 1.0 - risk_obj.lambda]
    end
    return params
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
    num_stages = get_num_stages(files)

    SAA = SDDPlab.Scenarios.generate_saa(scenarios, num_stages)

    return model, SAA
end
