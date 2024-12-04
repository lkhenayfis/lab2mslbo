using JuMP: JuMP
using SDDP: SDDP

function __default_ub_guess(files::Vector{SDDPlab.Inputs.InputModule}, stage::Integer)
    algo = SDDPlab.Inputs.get_algorithm(files)
    num_stages = SDDPlab.Inputs.get_number_of_stages(algo)

    buses = SDDPlab.System.get_buses_entities(SDDPlab.Inputs.get_system(files))
    scens = SDDPlab.Inputs.get_scenarios(files)

    guess = 0.0
    for stage_idx in range(stage, num_stages)
        for bus in buses
            guess += bus.deficit_cost * SDDPlab.Scenarios.get_load(bus.id, stage_idx, scens)
        end
    end

    return guess
end

function __default_α_guess(files::Vector{SDDPlab.Inputs.InputModule}, stage::Integer)
    algo = SDDPlab.Inputs.get_algorithm(files)
    num_stages = SDDPlab.Inputs.get_number_of_stages(algo)

    buses = SDDPlab.System.get_buses_entities(SDDPlab.Inputs.get_system(files))

    guess = 0.0
    for stage_idx in range(stage, num_stages)
        for bus in buses
            guess += bus.deficit_cost
        end
    end

    return guess
end

function __build_ub_model(
    files::Vector{SDDPlab.Inputs.InputModule}, optimizer
)::SDDP.PolicyGraph
    @info "Compiling upper model"
    sp_builder = __generate_subproblem_builder(files)
    algo = SDDPlab.Inputs.get_algorithm(files)
    num_stages = SDDPlab.Inputs.get_number_of_stages(algo)

    lip_builder(t) = __default_α_guess(files, t)
    ub_builder(t) = __default_ub_guess(files, t)
    ibf = InnerBellmanFunction(
        lip_builder; upper_bound = ub_builder, vertex_type = SDDP.SINGLE_CUT
    )

    pb_inner = SDDP.LinearPolicyGraph(
        sp_builder;
        stages = num_stages,
        sense = :Min,
        optimizer = optimizer,
        lower_bound = -Inf,
        upper_bound = Inf,
        bellman_function = ibf,
    )

    for (k, node) in pb_inner.nodes
        SDDP.set_objective(node)
    end

    return pb_inner
end

function __build_and_compute_ub_model(
    files::Vector{SDDPlab.Inputs.InputModule}, model::SDDP.PolicyGraph, optimizer
)
    @info "Evaluating upper policy"
    sp_builder = __generate_subproblem_builder(files)
    risk = SDDPlab.get_tasks(files)[2].risk_measure
    risk_measure = SDDPlab.Tasks.generate_risk_measure(risk)
    algo = SDDPlab.Inputs.get_algorithm(files)
    num_stages = SDDPlab.Inputs.get_number_of_stages(algo)

    lip_builder(t) = __default_α_guess(files, t)
    ub_builder(t) = __default_ub_guess(files, t)
    ibf = InnerBellmanFunction(
        lip_builder; upper_bound = ub_builder, vertex_type = SDDP.SINGLE_CUT
    )

    outer_model, upper_bound, upper_bound_time = build_compute_inner_dp(
        sp_builder,
        model;
        num_stages = num_stages,
        sense = :Min,
        optimizer = optimizer,
        lower_bound = 0.0,
        bellman_function = ibf,
        risk_measures = risk_measure,
    )

    return outer_model, upper_bound, upper_bound_time
end

function __generate_subproblem_builder(files::Vector{SDDPlab.Inputs.InputModule})::Function
    system = SDDPlab.Inputs.get_system(files)
    scenarios = SDDPlab.Inputs.get_scenarios(files)
    hydros_entities = SDDPlab.System.get_hydros_entities(system)
    num_hydros = length(hydros_entities)
    num_stages = SDDPlab.Inputs.get_number_of_stages(SDDPlab.Inputs.get_algorithm(files))

    SDDPlab.Scenarios.set_seed!(scenarios)

    SAA = SDDPlab.Scenarios.generate_saa(scenarios, num_stages)

    function fun_sp_build(m::JuMP.Model, node::Integer)
        SDDPlab.System.add_system_elements!(m, system)
        SDDPlab.Scenarios.add_uncertainties!(m, scenarios)

        # TODO - this will change once we have a proper load representation
        # as an stochastic process
        __add_load_balance!(m, files, node)

        Ω_node = vec(SAA[node])
        SDDP.parameterize(m, Ω_node) do ω
            JuMP.fix.(m[SDDPlab.Core.ω_INFLOW], ω)
            return nothing
        end

        SDDPlab.System.add_system_objective!(m, system)

        return nothing
    end

    return fun_sp_build
end

function __add_load_balance!(
    m::JuMP.Model, files::Vector{SDDPlab.Inputs.InputModule}, node::Integer
)
    system = SDDPlab.Inputs.get_system(files)
    hydros_entities = SDDPlab.System.get_hydros_entities(system)
    thermals_entities = SDDPlab.System.get_thermals_entities(system)
    lines_entities = SDDPlab.System.get_lines_entities(system)
    scenarios = SDDPlab.Inputs.get_scenarios(files)
    bus_ids = SDDPlab.System.get_ids(SDDPlab.System.get_buses(system))

    num_buses = length(bus_ids)
    num_lines = length(lines_entities)
    num_hydros = length(hydros_entities)
    num_thermals = length(thermals_entities)

    m[SDDPlab.Core.LOAD_BALANCE] = @constraint(
        m,
        [n = 1:num_buses],
        sum(
            m[SDDPlab.Core.HYDRO_GENERATION][j] for
            j in 1:num_hydros if hydros_entities[j].bus_id == bus_ids[n]
        ) +
        sum(
            m[SDDPlab.Core.THERMAL_GENERATION][j] for
            j in 1:num_thermals if thermals_entities[j].bus_id == bus_ids[n]
        ) +
        sum(
            m[SDDPlab.Core.DIRECT_EXCHANGE][j] - m[SDDPlab.Core.REVERSE_EXCHANGE][j] for
            j in 1:num_lines if lines_entities[j].target_bus_id == bus_ids[n]
        ) +
        sum(
            m[SDDPlab.Core.REVERSE_EXCHANGE][j] - m[SDDPlab.Core.DIRECT_EXCHANGE][j] for
            j in 1:num_lines if lines_entities[j].source_bus_id == bus_ids[n]
        ) +
        m[SDDPlab.Core.DEFICIT][bus_ids[n]] ==
            SDDPlab.Scenarios.get_load(bus_ids[n], node, scenarios)
    )
    return nothing
end

function __update_convergence_file(
    files::Vector{SDDPlab.Inputs.InputModule},
    upper_bound::Float64,
    upper_bound_time::Float64,
    e::CompositeException,
)
    tasks = SDDPlab.Inputs.Tasks.get_tasks(files)
    policy_index = findfirst(x -> isa(x, SDDPlab.Tasks.Policy), tasks)
    policy_output_path = tasks[policy_index].results.path
    policy_output_format = tasks[policy_index].results.format
    reader = SDDPlab.Inputs.Tasks.get_reader(policy_output_format)
    writer = SDDPlab.Inputs.Tasks.get_writer(policy_output_format)
    extension = SDDPlab.Inputs.Tasks.get_extension(policy_output_format)
    filename = policy_output_path * "/convergence" * extension
    convergence_df = reader(filename, e)
    # Adds upper_bound
    convergence_df[!, "upper_bound"] = fill(upper_bound, nrow(convergence_df))
    # Adds upper_bound_time
    rename!(convergence_df, "time" => "primal_time")
    convergence_df[!, "upper_bound_time"] = fill(0.0, nrow(convergence_df))
    convergence_df[nrow(convergence_df), "upper_bound_time"] = upper_bound_time
    convergence_df[!, "time"] =
        convergence_df[!, "primal_time"] + convergence_df[!, "upper_bound_time"]
    return writer(filename, convergence_df)
end
