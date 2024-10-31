
mutable struct SimpleLBOData
    x::Vector{VariableRef}
    c::Vector{Float64}
    A_eq::Matrix{Float64}
    B_eq::Matrix{Float64}
    T_eq::Matrix{Float64}
    d_eq::Vector{Vector{Float64}}
    A_ineq::Matrix{Float64}
    B_ineq::Matrix{Float64}
    T_ineq::Matrix{Float64}
    d_ineq::Vector{Vector{Float64}}
    Ux::Vector{Float64}
    Uy::Vector{Float64}
    lb::Float64
    ub::Float64
    α::Float64
    prob::Vector{Float64}

    SimpleLBOData() = new()
end

function SimpleLBOData(
    m::JuMP.Model, node_noises::Vector, lb::Float64, ub::Float64, α::Float64
)
    slbod = SimpleLBOData()

    md = ModelData(m)
    remove_dummy_variables!(md)
    fix_infinity_bounds!(md, 0.0, 1e7)
    force_inequalities_to_geq!(md)
    find_delete_ω!(md)

    fill_costs!(slbod, md)
    fill_equality_inequality_matrices!(slbod, md)
    fill_equality_inequality_rhs!(slbod, md, node_noises)
    fill_bounds!(slbod, md)

    slbod.x = md.x
    slbod.lb = lb
    slbod.ub = ub
    slbod.α = α

    return slbod
end

# HELPERS ------------------------------------------------------------------------------------------

function fill_costs!(slbod::SimpleLBOData, md::ModelData)
    state_t, state_t_1, control = get_state_control_indexes(md)
    return slbod.c = md.c[control]
end

function fill_equality_inequality_matrices!(slbod::SimpleLBOData, md::ModelData)
    eq_ind, ineq_ind = get_equality_inequality_indexes(md)
    A, B, T = split_constraint_matrices(md)

    slbod.A_eq = A[eq_ind, :]
    slbod.B_eq = B[eq_ind, :]
    slbod.T_eq = T[eq_ind, :]

    slbod.A_ineq = A[ineq_ind, :]
    slbod.B_ineq = B[ineq_ind, :]
    return slbod.T_ineq = T[ineq_ind, :]
end

function fill_equality_inequality_rhs!(
    slbod::SimpleLBOData,
    md::ModelData,
    noises::Vector;
    prob::Vector{Float64} = ones(size(noises)) / length(noises),
)
    eq_ind, ineq_ind = get_equality_inequality_indexes(md)

    M = length(noises)
    ds = [copy(md.b_lower) for i in range(1, M)]
    for (i, d) in enumerate(ds)
        d[md.sp_constraints] .= noises[i]
    end

    slbod.d_eq = [d[eq_ind] for d in ds]
    slbod.d_ineq = [d[ineq_ind] for d in ds]
    return slbod.prob = prob
end

function fill_bounds!(slbod::SimpleLBOData, md::ModelData)
    state_t, state_t_1, control = get_state_control_indexes(md)

    lb_states = md.x_lower[state_t]
    ub_states = md.x_upper[state_t]

    lb_control = md.x_lower[control]
    ub_control = md.x_upper[control]

    slbod.Ux = ub_states
    return slbod.Uy = ub_control
end