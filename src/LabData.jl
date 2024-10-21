mutable struct LabData
    state0::Vector{Float64}
    seed::Int64
    num_stages::Int64
    num_iterations::Int64
    risk_alpha::Float64
    risk_lambda::Float64
end
