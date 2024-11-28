using SDDPlab: SDDPlab
using lab2mslbo: lab2mslbo

deck_dir = "./data-1dtoy/"

e = CompositeException()
curdir = pwd()
cd(deck_dir)

# Runs policy evaluation
entrypoint = SDDPlab.Inputs.Entrypoint("main.jsonc", e)
artifacts = SDDPlab.__run_tasks!(entrypoint, e)
SDDPlab.__log_errors(e)

# Build and eval ub model
policy_index = findfirst(x -> isa(x, SDDPlab.Tasks.PolicyArtifact), artifacts)
policy = artifacts[policy_index].policy

# Transforms to vertex policy graph
inner_policy, upper_bound, upper_bound_time = lab2mslbo.__build_and_compute_ub_model(
    entrypoint.inputs.files, policy
)

lab2mslbo.__update_convergence_file(
    entrypoint.inputs.files, upper_bound, upper_bound_time, e
)

# Generates fake policy artifact and artifact vector
task_definitions = SDDPlab.get_tasks(entrypoint.inputs.files)
policy_task_index = findfirst(x -> isa(x, SDDPlab.Tasks.Policy), task_definitions)
policy_task_definition = task_definitions[policy_task_index]

artifacts = Vector{SDDPlab.Tasks.TaskArtifact}([
    SDDPlab.Tasks.InputsArtifact(entrypoint.inputs.path, entrypoint.inputs.files),
    SDDPlab.Tasks.PolicyArtifact(
        policy_task_definition, inner_policy, entrypoint.inputs.files
    ),
])

# Runs simulation again
simulation_task_index = findfirst(x -> isa(x, SDDPlab.Tasks.Simulation), task_definitions)
simulation_task_definition = task_definitions[simulation_task_index]
a = SDDPlab.Tasks.run_task(simulation_task_definition, artifacts, e)
SDDPlab.__save_results(a)