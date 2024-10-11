using DualSDDP
using Random: seed!
import GLPK

solver = GLPK.Optimizer

include("src/lab2mslbo.jl")
M = lab2mslbo.build_mslbo("data-1dtoy/");

state0 = [83.222]
risk      = mk_primal_avar(.4)
risk_dual = mk_copersp_avar(.4)


# Need M::MSLBO, nstages::Int, risk, risk_dual, state0::Vector{Float64},
#      niters::Int, solver
# Solution algorithms
# Pure primal
seed!(2)
primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(M, 4, risk, solver, state0, 40; verbose=true)

# Pure dual
seed!(1)
dual_pb, dual_ubs, dual_times = dualsolve(M, 4, risk_dual, solver, state0, 40; verbose=true)

# Recursive upper bounds over primal trajectories
rec_ubs, rec_times = primalub(M, 4, risk, solver, primal_trajs, 10:10:40; verbose=true)

# Primal with outer and inner bounds
seed!(1)
io_pb, io_lbs, io_ubs, io_times = problem_child_solve(M, 4, risk, solver, state0, 40; verbose=true)

