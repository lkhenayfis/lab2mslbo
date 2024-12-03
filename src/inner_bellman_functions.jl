#  Copyright (c) 2024 Bernardo Freitas Paulo da Costa
#
#  This code is heavily based on "bellman_functions.jl" from SDDP.jl
#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

using JuMP: JuMP
using SDDP: SDDP
using JuMP: @variable, @constraint, @expression, set_normalized_coefficient

mutable struct Vertex
    value::Float64
    state::Dict{Symbol,Float64}
    obj_y::Union{Nothing,NTuple{N,Float64} where {N}}
    belief_y::Union{Nothing,Dict{T,Float64} where {T}}
    non_dominated_count::Int
    variable_ref::Union{Nothing,JuMP.VariableRef}
end

mutable struct InnerConvexApproximation
    theta::JuMP.VariableRef
    states::Dict{Symbol,JuMP.VariableRef}
    objective_states::Union{Nothing,NTuple{N,JuMP.VariableRef} where {N}}
    belief_states::Union{Nothing,Dict{T,JuMP.VariableRef} where {T}}
    # Storage for cut selection
    vertices::Vector{Vertex}
    sampled_states::Vector{SDDP.SampledState}
    vertices_to_be_deleted::Vector{Vertex}
    deletion_minimum::Int
    Lipschitz_constant::Float64

    function InnerConvexApproximation(
        theta::JuMP.VariableRef,
        states::Dict{Symbol,JuMP.VariableRef},
        objective_states,
        belief_states,
        deletion_minimum::Int,
        Lipschitz_constant::Float64,
    )
        return new(
            theta,
            states,
            objective_states,
            belief_states,
            Vertex[],
            SDDP.SampledState[],
            Vertex[],
            deletion_minimum,
            Lipschitz_constant,
        )
    end
end

function _add_vertex(
    V::InnerConvexApproximation,
    θᵏ::Float64,
    xᵏ::Dict{Symbol,Float64},
    obj_y::Union{Nothing,NTuple{N,Float64}},
    belief_y::Union{Nothing,Dict{T,Float64}};
    cut_selection::Bool = true,
) where {N,T}
    vertex = Vertex(θᵏ, xᵏ, obj_y, belief_y, 1, nothing)
    _add_vertex_var_to_model(V, vertex)
    push!(V.vertices, vertex)
    if cut_selection
        # TODO(bfpc): implement cut selection
        nothing # _cut_selection_update(V, vertex, xᵏ)
    end
    return nothing
end

# TODO(bfpc): implement Objective and Belief states
function _add_vertex_var_to_model(V::InnerConvexApproximation, vertex::Vertex)
    model = JuMP.owner_model(V.theta)
    yᵀμ = JuMP.AffExpr(0.0)
    if V.objective_states !== nothing
        error("Objective states not yet implemented.")
        for (y, μ) in zip(cut.obj_y, V.objective_states)
            JuMP.add_to_expression!(yᵀμ, y, μ)
        end
    end
    if V.belief_states !== nothing
        error("Belief states not yet implemented.")
        for (k, μ) in V.belief_states
            JuMP.add_to_expression!(yᵀμ, cut.belief_y[k], μ)
        end
    end

    # Add a new variable to the convex combination constraints
    σk = @variable(model, lower_bound = 0.0, upper_bound = 1.0)
    set_normalized_coefficient(model[:σ_cc], σk, 1)

    xk = vertex.state
    for key in keys(xk)
        set_normalized_coefficient(model[:x_cc][key], σk, xk[key])
    end

    vk = vertex.value
    set_normalized_coefficient(model[:theta_cc], σk, vk)

    vertex.variable_ref = σk
    return nothing
end

"""
Exact vertex selection

Test if (x_j, v_j) is covered by (cc x_i, cc v_i) + (0, t)

We could add the Lipschitz cone, generated by [+ (±e_j, Lip)],
but in practice this should almost never happen.
"""
function _vertex_selection(V::InnerConvexApproximation, solver)
    nvertices = length(V.vertices)
    state_keys = keys(V.states)

    m = JuMP.Model(solver)
    JuMP.set_silent(m)
    @variable(m, 0 <= σ[1:nvertices] <= 1)
    @variable(m, t)
    @objective(m, Max, t)
    @constraint(
        m,
        cc_x[k in state_keys],
        sum(σ[i] * V.vertices[i].state[k] for i in 1:nvertices) .== 0
    )
    @constraint(m, cc_t, sum(σ[i] * V.vertices[i].value for i in 1:nvertices) + t == 0)
    @constraint(m, sum(σ[i] for i in 1:nvertices) == 1)

    remove_vars = VariableRef[]
    for j in 1:nvertices
        v = V.vertices[j]
        for k in state_keys
            JuMP.set_normalized_rhs(cc_x[k], v.state[k])
        end
        JuMP.set_normalized_rhs(cc_t, v.value)

        optimize!(m)
        margin = objective_value(m)
        if margin > 1e-6
            # println("Vertex $(v) ($j-th) is $(margin) redundant")
            # println(value.(σ))
            push!(remove_vars, v.variable_ref)
            # Don't use this vertex anymore for vertex selection
            delete(m, σ[j])
        end
    end
    # Remove all vertices from subproblem
    if length(remove_vars) > 0
        sp = owner_model(V.theta)
        delete(sp, remove_vars)
        println("Selection removed $(length(remove_vars)) vertices")
    end

    return m
end

"""
    InnerBellmanFunction

An inner representation of the value function, currently only using

 1. x - convex "resource" states

This function penalizes the convex combination constraint to define a valid
approximation in the entire state space.

In addition, we have three types of vertices (x, v) (in "cut" terminology):

 1. Single-cuts (also called "average" cuts in the literature), which involve
    the risk-adjusted expectation of the cost-to-go.
 2. Multi-cuts, which use a different cost-to-go term for each realization w.
 3. Risk-cuts, which correspond to the facets of the dual interpretation of a
    convex risk measure.

Therefore, InnerBellmanFunction returns a JuMP model of the following form (when minimizing):

```
V(x, b, y) =
    min: θ
    s.t. # "Single" / "Average" cuts
         θ ≥ Σ α(j) * v(j) + Lip |delta|
         x = Σ α(j) * x(j) + delta
         1 = Σ α(j)
         0 ≤ α(j) ≤ 1,                                ∀ j ∈ J
         # "Multi" cuts
         φ(w) ≥ Σ α(k, w) * v(k, w) + Lip |delta(w)|, ∀ w ∈ Ω
         x = Σ α(k, w) * x(k, w) + delta(w),          ∀ w ∈ Ω
         1 = Σ α(k, w),                               ∀ w ∈ Ω
         0 ≤ α(k, w) ≤ 1,                             ∀ k ∈ K, w ∈ Ω
         # "Risk-set" cuts
         θ ≥ Σ{p(k, w) * φ(w)}_w - μᵀb(k) - νᵀy(k),   ∀ k ∈ K
```
"""
mutable struct InnerBellmanFunction
    cut_type::SDDP.CutType
    global_theta::InnerConvexApproximation
    local_thetas::Vector{InnerConvexApproximation}
    # Cuts defining the dual representation of the risk measure.
    risk_set_cuts::Set{Vector{Float64}}
    Lipschitz_constant::Float64
end

"""
    InnerBellmanFunction(Lipschitz_constant;
        lower_bound = -Inf,
        upper_bound = Inf,
        deletion_minimum::Int = 1,
        vertex_type::SDDP.CutType = SDDP.MULTI_CUT,
    )
"""
function InnerBellmanFunction(
    Lipschitz_constant;
    lower_bound = -Inf,
    upper_bound = Inf,
    deletion_minimum::Int = 1,
    vertex_type::SDDP.CutType = SDDP.MULTI_CUT,
)
    return SDDP.InstanceFactory{InnerBellmanFunction}(;
        Lipschitz_constant = Lipschitz_constant,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        deletion_minimum = deletion_minimum,
        vertex_type = vertex_type,
    )
end

function SDDP.bellman_term(bellman_function::InnerBellmanFunction)
    return bellman_function.global_theta.theta
end

# to_node is an internal helper function so users can pass arguments like:
#   Lipschitz_constant = 2000.0,
#   Lipschitz_constant = Dict(1=>2000.0, 2=>1000.0)
#   Lipschitz_constant = (node_index) -> 1000.0 * (nStages - node_index)
# See [`to_nodal_form`](@ref)
function to_node(node::SDDP.Node{T}, element) where {T}
    return element
end

function to_node(node::SDDP.Node{T}, builder::Function) where {T}
    return builder(node.index)
end

function to_node(node::SDDP.Node{T}, dict::Dict{T,V}) where {T,V}
    return dict[node.index]
end

function SDDP.initialize_bellman_function(
    factory::SDDP.InstanceFactory{InnerBellmanFunction},
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
) where {T}
    lower_bound, upper_bound, deletion_minimum, vertex_type, Lipschitz_constant = -Inf,
    Inf, 0, SDDP.SINGLE_CUT,
    0.0
    if length(factory.args) > 0
        error("Positional arguments $(factory.args) ignored in InnerBellmanFunction.")
    end
    # TODO(bfpc): generalize to dictionaries for lb/ub/Lip, taking node as keys
    for (kw, value) in factory.kwargs
        value = to_node(node, value)
        if kw == :lower_bound
            lower_bound = value
        elseif kw == :upper_bound
            upper_bound = value
        elseif kw == :deletion_minimum
            deletion_minimum = value
        elseif kw == :vertex_type
            vertex_type = value
        elseif kw == :Lipschitz_constant
            Lipschitz_constant = value
        else
            error("Keyword $(kw) not recognised as argument to InnerBellmanFunction.")
        end
    end
    if lower_bound == -Inf && upper_bound == Inf
        error("You must specify a finite bound on the inner cost-to-go term.")
    end
    if length(node.children) == 0
        lower_bound = upper_bound = 0.0
    end

    # Add epigraph variable for inner approximation
    sp = node.subproblem
    Θᴳ = @variable(sp)
    # Initialize bounds for the objective states. If objective_state==nothing,
    # this check will be skipped by dispatch.
    SDDP._add_initial_bounds(node.objective_state, Θᴳ)
    x′ = Dict(key => var.out for (key, var) in node.states)
    obj_μ = node.objective_state !== nothing ? node.objective_state.μ : nothing
    belief_μ = node.belief_state !== nothing ? node.belief_state.μ : nothing

    # Model convex combination + Lipschitz penalty
    if length(node.children) != 0
        @variable(sp, δ[keys(x′)])
        @variable(sp, δ_abs[keys(x′)] >= 0)
        @constraint(sp, [k in keys(x′)], δ_abs[k] >= δ[k])
        @constraint(sp, [k in keys(x′)], δ_abs[k] >= -δ[k])

        @variable(sp, 0 <= σ0 <= 1)
        if JuMP.objective_sense(sp) == JuMP.MOI.MIN_SENSE
            @constraint(
                sp, theta_cc, Lipschitz_constant * sum(δ_abs) + σ0 * upper_bound <= Θᴳ
            )
        else
            @constraint(
                sp, theta_cc, Lipschitz_constant * sum(δ_abs) + σ0 * lower_bound >= Θᴳ
            )
        end
        @constraint(sp, x_cc[k in keys(x′)], δ[k] == x′[k])
        @constraint(sp, σ_cc, σ0 == 1)

        @expression(sp, VERTEX_COVERAGE_DISTANCE, sum(δ_abs))

    else
        JuMP.fix(Θᴳ, 0.0)
        @expression(sp, VERTEX_COVERAGE_DISTANCE, 0.0)
    end

    return InnerBellmanFunction(
        vertex_type,
        InnerConvexApproximation(
            Θᴳ, x′, obj_μ, belief_μ, deletion_minimum, Lipschitz_constant
        ),
        InnerConvexApproximation[],
        Set{Vector{Float64}}(),
        Lipschitz_constant,
    )
end

function refine_inner_bellman_function(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    bellman_function::InnerBellmanFunction,
    risk_measure::SDDP.AbstractRiskMeasure,
    outgoing_state::Dict{Symbol,Float64},
    dual_variables::Vector{Dict{Symbol,Float64}},
    noise_supports::Vector,
    nominal_probability::Vector{Float64},
    objective_realizations::Vector{Float64},
) where {T}
    lock(node.lock)
    try
        return _refine_inner_bellman_function_no_lock(
            model,
            node,
            bellman_function,
            risk_measure,
            outgoing_state,
            dual_variables,
            noise_supports,
            nominal_probability,
            objective_realizations,
        )
    finally
        unlock(node.lock)
    end
end

function _refine_inner_bellman_function_no_lock(
    model::SDDP.PolicyGraph{T},
    node::SDDP.Node{T},
    bellman_function::InnerBellmanFunction,
    risk_measure::SDDP.AbstractRiskMeasure,
    outgoing_state::Dict{Symbol,Float64},
    dual_variables::Vector{Dict{Symbol,Float64}},
    noise_supports::Vector,
    nominal_probability::Vector{Float64},
    objective_realizations::Vector{Float64},
) where {T}
    # Sanity checks.
    @assert length(dual_variables) ==
        length(noise_supports) ==
        length(nominal_probability) ==
        length(objective_realizations)
    # Preliminaries that are common to all cut types.
    risk_adjusted_probability = similar(nominal_probability)
    offset = SDDP.adjust_probability(
        risk_measure,
        risk_adjusted_probability,
        nominal_probability,
        noise_supports,
        objective_realizations,
        model.objective_sense == JuMP.MOI.MIN_SENSE,
    )
    # The meat of the function.
    if bellman_function.cut_type == SDDP.SINGLE_CUT
        return _add_average_vertex(
            node,
            outgoing_state,
            risk_adjusted_probability,
            objective_realizations,
            dual_variables,
            offset,
        )
    else  # TODO(bfpc): Add a multi-cut
        error("Multi-vertices not yet implemented.")
        @assert bellman_function.cut_type == SDDP.MULTI_CUT
        SDDP._add_locals_if_necessary(node, bellman_function, length(dual_variables))
        return _add_multi_vertex(
            node,
            outgoing_state,
            risk_adjusted_probability,
            objective_realizations,
            dual_variables,
            offset,
        )
    end
end

function _add_average_vertex(
    node::SDDP.Node,
    outgoing_state::Dict{Symbol,Float64},
    risk_adjusted_probability::Vector{Float64},
    objective_realizations::Vector{Float64},
    dual_variables::Vector{Dict{Symbol,Float64}},
    offset::Float64,
)
    N = length(risk_adjusted_probability)
    @assert N == length(objective_realizations) == length(dual_variables)
    # Calculate the expected intercept with respect to the
    # risk-adjusted probability distribution.
    θᵏ = offset
    for i in 1:length(objective_realizations)
        p = risk_adjusted_probability[i]
        θᵏ += p * objective_realizations[i]
    end
    # Now add the average-vertex to the subproblem. We include the objective-state
    # component μᵀy and the belief state (if it exists).
    obj_y = node.objective_state === nothing ? nothing : node.objective_state.state
    belief_y = node.belief_state === nothing ? nothing : node.belief_state.belief
    _add_vertex(node.bellman_function.global_theta, θᵏ, outgoing_state, obj_y, belief_y)
    return (theta = θᵏ, x = outgoing_state, obj_y = obj_y, belief_y = belief_y)
end

function model_type(model::SDDP.PolicyGraph{T}) where {T}
    return T
end

function build_compute_inner_dp(
    build::Function,
    pb::SDDP.PolicyGraph;
    num_stages::Int,
    sense::Symbol = :Min,
    optimizer,
    lower_bound::Float64 = -Inf,
    upper_bound::Float64 = Inf,
    bellman_function,
    risk_measures::SDDP.AbstractRiskMeasure,
    print_level::Int = 1,
)
    pb_inner = SDDP.LinearPolicyGraph(
        build;
        stages = num_stages,
        sense,
        optimizer,
        lower_bound,
        upper_bound,
        bellman_function,
    )

    for (k, node) in pb_inner.nodes
        SDDP.set_objective(node)
    end

    opts = SDDP.Options(pb_inner, pb.initial_root_state; risk_measures)
    T = model_type(pb_inner)
    total_dt = 0.0
    for node_index in sort(collect(keys(pb.nodes)); rev = true)[2:end]
        dt = @elapsed begin
            node = pb_inner[node_index]
            fw_samples = pb[node_index].bellman_function.global_theta.sampled_states
            for sampled_state in fw_samples
                outgoing_state = sampled_state.state
                items = SDDP.BackwardPassItems(T, SDDP.Noise)
                SDDP.solve_all_children(
                    pb_inner,
                    node,
                    items,
                    1.0,
                    nothing, # belief_state
                    nothing, # objective_state
                    outgoing_state,
                    opts.backward_sampling_scheme,
                    Tuple{T,Any}[], # scenario_path
                    opts.duality_handler,
                    opts,
                )
                refine_inner_bellman_function(
                    pb_inner,
                    node,
                    node.bellman_function,
                    opts.risk_measures[node_index],
                    outgoing_state,
                    items.duals,
                    items.supports,
                    items.probability,
                    items.objectives,
                )
            end
        end
        dt_vs = @elapsed _vertex_selection(node.bellman_function.global_theta, optimizer)
        if print_level > 0
            println(
                "Node: $(node_index) - elapsed time: ",
                round(dt; digits = 2),
                " plus ",
                round(dt_vs; digits = 2),
                " for vertex selection.",
            )
        end
        total_dt += dt + dt_vs
    end

    ub = SDDP.calculate_bound(pb_inner; risk_measure = risk_measures)
    if print_level > 0
        println("First-stage upper bound: ", ub)
        println("Total time for upper bound: ", total_dt)
    end
    return pb_inner, ub, total_dt
end

function _get_vertex_info(
    json_vertex, dualcuts; vertex_name_parser::Function = _vertex_name_parser
)
    if dualcuts
        value = -json_vertex["intercept"]
        state = Dict(
            Symbol(vertex_name_parser(k)) => v for (k, v) in json_vertex["coefficients"]
        )
    else
        value = json_vertex["value"]
        state = Dict(Symbol(vertex_name_parser(k)) => v for (k, v) in json_vertex["state"])
    end
    obj_y = nothing
    belief_y = nothing
    return value, state, obj_y, belief_y
end

function _vertex_name_parser(vertex_name::String)::String
    return vertex_name
end

"""
    read_vertices_from_file(
        model::PolicyGraph{T},
        filename::String;
        kwargs...,
    ) where {T}

Read vertices from `filename` into `model`.

Since `T` can be an arbitrary Julia type, the conversion to JSON is lossy. When
reading, `read_cuts_from_file` only supports `T=Int`, `T=NTuple{N, Int}`, and
`T=Symbol`. If you have manually created a policy graph with a different node
type `T`, provide a function `node_name_parser` with the signature

## Keyword arguments

  - `dualcuts::Bool` transform dual cuts into vertices.

  - `node_name_parser(T, name::String)::T where {T}` that returns the name of each
    node given the string name `name`.
    If `node_name_parser` returns `nothing`, those cuts are skipped.
  - `vertex_selection::Bool` run or not the exact vertex selection algorithm when
    adding vertices to the model.
"""
function read_vertices_from_file(
    model::SDDP.PolicyGraph{T},
    filename::String;
    dualcuts::Bool = false,
    node_name_parser::Function = (::Type{Int64}, x::String) -> parse(Int64, x),
    vertex_name_parser::Function = _vertex_name_parser,
    vertex_selection::Bool = true,
) where {T}
    vertices = JSON.parsefile(filename; use_mmap = false)
    for node_info in vertices
        node_name = node_name_parser(T, node_info["node"])::Union{Nothing,T}
        if node_name === nothing
            continue  # Skip reading these vertices
        end
        node = model[node_name]
        bf = node.bellman_function
        # Loop through and add the vertices.
        list = (dualcuts ? node_info["single_cuts"] : node_info["vertices"])
        for json_vertex in list
            value, state, obj_y, belief_y = _get_vertex_info(
                json_vertex, dualcuts; vertex_name_parser = vertex_name_parser
            )
            _add_vertex(bf.global_theta, value, state, obj_y, belief_y)
        end
        if vertex_selection
            _vertex_selection(bf.global_theta, optimizer)
        end

        # Loop through and add the multi-cuts. There are two parts:
        #  (i) the cuts w.r.t. the state variable x
        # (ii) the cuts that define the risk set
        # There is one additional complication: if these cuts are being read
        # into a new model, the local theta variables may not exist yet.
        if length(node_info["risk_set_cuts"]) > 0
            _add_locals_if_necessary(node, bf, length(first(node_info["risk_set_cuts"])))
        end
        list = (dualcuts ? node_info["multi_cuts"] : node_info["multi_vertices"])
        for json_vertex in list
            @error "Multi-vertices not yet implemented."
            value, state, obj_y, belief_y = _get_vertex_info(json_vertex, dualcuts)
            _add_vertex(
                bf.local_thetas[json_vertex["realization"]], value, state, obj_y, belief_y
            )
        end
        if vertex_selection && length(list) > 0
            for local_vf in bf.local_thetas
                _vertex_selection(local_vf, optimizer)
            end
        end
        # Here is part (ii): adding the constraints that define the risk-set
        # representation of the risk measure.
        for json_cut in node_info["risk_set_cuts"]
            expr = @expression(
                node.subproblem,
                bf.global_theta.theta -
                    sum(p * V.theta for (p, V) in zip(json_cut, bf.local_thetas))
            )
            if JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE
                @constraint(node.subproblem, expr >= 0)
            else
                @constraint(node.subproblem, expr <= 0)
            end
        end
    end
    return nothing
end
