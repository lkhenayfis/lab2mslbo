using DualSDDP
using Random: seed!
import GLPK

solver = GLPK.Optimizer

include("src/lab2mslbo.jl")
M, data = lab2mslbo.build_mslbo("data-1dtoy/");

state0 = data.state0
risk = mk_primal_avar(data.risk_alpha; beta=data.risk_lambda)
risk_dual = mk_copersp_avar(data.risk_alpha; beta=data.risk_lambda)
nstages = data.num_stages
niters = data.num_iterations
ub_iters = Int64.(2 .^ (0:1:floor(log2(niters))))

# O que fazer com os seeds?
# Need solver

# Solution algorithms

# Pure primal
seed!(2)
primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(M, nstages, risk, solver, state0, niters; verbose=true);

# Pure dual
seed!(1)
dual_pb, dual_ubs, dual_times = dualsolve(M, nstages, risk_dual, solver, state0, niters; verbose=true);

# Recursive upper bounds over primal trajectories
rec_ubs, rec_times = primalub(M, nstages, risk, solver, primal_trajs, ub_iters; verbose=true);

# Primal with outer and inner bounds
seed!(1)
io_pb, io_lbs, io_ubs, io_times = problem_child_solve(M, nstages, risk, solver, state0, niters; verbose=true);

