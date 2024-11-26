using DualSDDP
using SDDPlab: SDDPlab
using Random: seed!
using lab2mslbo: lab2mslbo
using DataFrames

deck_dir = "./data-1dtoy/"
curdir = pwd()

## ---------- InnerBellmanFunction calls ------------

cd(deck_dir)
e = CompositeException()

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

# Generates fake policy artifact
# Index 1 has a EchoArtifact, Index 2 a PolicyArtifact...
artifacts[2] = SDDPlab.Tasks.PolicyArtifact(
    artifacts[policy_index].task, inner_policy, entrypoint.inputs.files
)

# Runs simulation again
simulation_index = findfirst(x -> isa(x, SDDPlab.Tasks.SimulationArtifact), artifacts)
simulation = artifacts[simulation_index].task
a = SDDPlab.Tasks.run_task(simulation, artifacts, e)
SDDPlab.__save_results(a)

## ---------- DualSDDP calls ------------

cd(curdir)

M, data = lab2mslbo.build_mslbo(deck_dir)

mkpath(deck_dir * data.output_path)

risk = mk_primal_avar(data.risk_alpha; beta = data.risk_lambda)
risk_dual = mk_copersp_avar(data.risk_alpha; beta = data.risk_lambda)
ub_iters = Int64.(2 .^ (4:1:floor(log2(data.num_iterations))))

## Solution algorithms

# Pure primal
seed!(data.seed)
primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(
    M, data.num_stages, risk, data.solver, data.state0, data.num_iterations; verbose = true
);

# Pure dual
seed!(data.seed)
dual_pb, dual_ubs, dual_times = dualsolve(
    M,
    data.num_stages,
    risk_dual,
    data.solver,
    data.state0,
    data.num_iterations;
    verbose = true,
);

# Recursive upper bounds over primal trajectories
rec_ubs, rec_times = primalub(
    M, data.num_stages, risk, data.solver, primal_trajs, ub_iters; verbose = true
);

# Primal with outer and inner bounds
seed!(data.seed)
io_pb, io_lbs, io_ubs, io_times = problem_child_solve(
    M, data.num_stages, risk, data.solver, data.state0, data.num_iterations; verbose = true
);

## Exporting outputs

lab2mslbo.export_primal_with_ub_convergence(
    data.num_iterations,
    primal_lbs,
    primal_times,
    rec_ubs,
    rec_times,
    data.writer,
    data.extension;
    output_path_without_extension = deck_dir * data.output_path * "/convergence_philpott",
)

DualSDDP.write_cuts_to_file(primal_pb, deck_dir * data.output_path * "/primal_cuts.json")

lab2mslbo.export_primal_with_dual_ub_convergence(
    data.num_iterations,
    primal_lbs,
    primal_times,
    dual_ubs,
    dual_times,
    data.writer,
    data.extension;
    output_path_without_extension = deck_dir * data.output_path * "/convergence_dual",
)

DualSDDP.write_cuts_to_file(dual_pb, deck_dir * data.output_path * "/dual_cuts.json")

lab2mslbo.export_problem_child_convergence(
    data.num_iterations,
    io_lbs,
    io_ubs,
    io_times,
    data.writer,
    data.extension;
    output_path_without_extension = deck_dir * data.output_path * "/convergence_reagan",
)

DualSDDP.write_policy_to_file(io_pb, deck_dir * data.output_path * "/reagan_policy.json")