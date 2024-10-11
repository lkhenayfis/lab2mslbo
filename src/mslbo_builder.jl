
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
    seed::Integer = 1234)::Tuple{DualSDDP.MSLBO, Vector{Float64}}

    files = read_lab_inputs(data_dir)
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

    return DualSDDP.build(lbos), initial_states
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

function get_initial_states(files::Vector{SDDPlab.Inputs.InputModule})
    hydros = SDDPlab.System.get_hydros_entities(SDDPlab.Inputs.get_system(files))
    initial_states = [h.initial_storage for h in hydros]

    return initial_states
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
