
mutable struct SimpleLBOData
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

function SimpleLBOData(m::JuMP.Model, node_noises::Vector)

    slbod = SimpleLBOData()

    md = ModelData(m)
    remove_dummy_variables!(md)
    fix_infinity_bounds!(md)

    # preencher custos
        # fill_costs!(slbod, md)
    # preencher matrizes A B T de igualdade e desigualdade
    # preencher bounds de estados e controle

end