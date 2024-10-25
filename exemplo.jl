using DualSDDP
using SDDPlab: SDDPlab
using Random: seed!
using lab2mslbo: lab2mslbo
using DataFrames


deck_dir = "./data-1dtoy/"

curdir = pwd()
cd(deck_dir)
SDDPlab.main()
cd(curdir)

M, data, solver, writer, extension = lab2mslbo.build_mslbo(deck_dir);

output_path = data.output_path
state0 = data.state0
seed = data.seed
risk = mk_primal_avar(data.risk_alpha; beta=data.risk_lambda)
risk_dual = mk_copersp_avar(data.risk_alpha; beta=data.risk_lambda)
nstages = data.num_stages
niters = data.num_iterations
ub_iters = Int64.(2 .^ (0:1:floor(log2(niters))))


# Solution algorithms

# Pure primal
seed!(seed)
primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(M, nstages, risk, solver, state0, niters; verbose=true);

# Pure dual
seed!(seed)
dual_pb, dual_ubs, dual_times = dualsolve(M, nstages, risk_dual, solver, state0, niters; verbose=true);

# Recursive upper bounds over primal trajectories
rec_ubs, rec_times = primalub(M, nstages, risk, solver, primal_trajs, ub_iters; verbose=true);

# Primal with outer and inner bounds
seed!(seed)
io_pb, io_lbs, io_ubs, io_times = problem_child_solve(M, nstages, risk, solver, state0, niters; verbose=true);


# Convergence columns - sddp-lab
# iteration | lower_bound | simulation | upper_bound | time

# Philpott
# iteration | primal_lbs  |            |   rec_ubs   | primal_times | rec_times | time

dense_ubs = fill(NaN, niters)
for ub_pair in rec_ubs
    dense_ubs[ub_pair[1]] = ub_pair[2]
end
dense_rec_times = fill(NaN, niters)
for (ub_pair, time) in zip(rec_ubs, rec_times)
    dense_rec_times[ub_pair[1]] = time
end

philpott_df = DataFrame(iteration=1:niters,
    lower_bound=primal_lbs,
    simulation=fill(NaN, niters),
    upper_bound=dense_ubs,
    primal_time=primal_times,
    upper_bound_time=dense_rec_times,
    time=primal_times + replace(dense_rec_times, NaN => 0.0))


output_dir_path = deck_dir * output_path
mkpath(output_dir_path)
writer(output_dir_path * "/convergence_philpott" * extension, philpott_df)

# DualSDDP
# iteration | primal_lbs  |            |   dual_ubs  | primal_times | dual_times | time

dual_df = DataFrame(iteration=1:niters,
    lower_bound=primal_lbs,
    simulation=fill(NaN, niters),
    upper_bound=dual_ubs,
    primal_time=primal_times,
    upper_bound_time=dual_times,
    time=primal_times + dual_times)

output_dir_path = deck_dir * output_path
mkpath(output_dir_path)
writer(output_dir_path * "/convergence_dual" * extension, dual_df)


# Reagan (Baucke)
# iteration | io_lbs      |            |   io_ubs    |   io_times   | time

reagan_df = DataFrame(iteration=1:niters,
    lower_bound=io_lbs,
    simulation=fill(NaN, niters),
    upper_bound=io_ubs,
    time=io_times)

output_dir_path = deck_dir * output_path
mkpath(output_dir_path)
writer(output_dir_path * "/convergence_reagan" * extension, reagan_df)

