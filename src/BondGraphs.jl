module BondGraphs

using DifferentialEquations
using ModelingToolkit
using SymbolicUtils
using Symbolics
using LinearAlgebra
using Graphs
using MetaGraphs
using FileIO
using Documenter
## TODO: Develop better model parameter handling avoid structs of structs....
## TODO: Better interface with LightGraphs and Flux
export BondGraph, BioBondGraph,
    add_Bond!, add_R!, add_C!, add_I!, add_M!,
    add_Se!, add_Sf!,
    add_TF!, add_GY!,
    add_MTF!, add_MGY!,
    add_1J!, add_0J!,
    add_C_multiport!, add_I_multiport!,
    generate_model!, generate_model,
    state_space,
    get_graph, savebondgraph, loadbondgraph,
    remove_algebraic

## Function to create a generic Model
abstract type AbstractBondGraph end

mutable struct BondGraph <: AbstractBondGraph
    graph::MetaDiGraph
    model::ODESystem
end

mutable struct BioBondGraph <: AbstractBondGraph
    graph::MetaDiGraph
    model::ODESystem
    R
    T
end
"""

Get the ODE System Corresponding to the Specific Element

"""
Base.getindex(BG::AbstractBondGraph, node::Symbol) = get_prop(BG.graph, BG.graph[node, :name], :sys)

"""

Get the system corresponding to Node - `node`

"""
Base.getindex(BG::AbstractBondGraph, node::Int) = get_prop(BG.graph, node, :sys)


"""

Create an empty BondGraph to be populated during the analysis

"""
function BondGraph(independent_variable::Num)
    mg = MetaDiGraph()
    set_indexing_prop!(mg, :name)
    sys = ODESystem(Equation[], independent_variable, name = :model)
    return BondGraph(mg, sys)
end

function BioBondGraph(independent_variable::Num; R = 1.0, T = 1.0)
    mg = MetaDiGraph()
    set_indexing_prop!(mg, :name)
    sys = ODESystem(Equation[], independent_variable, name = :model)
    return BioBondGraph(mg, sys, R, T)
end

"""

Create a BondGraph provided a directed metagraph with node `:name` and `:type` defined for each node

"""
function BondGraph(mg::AbstractMetaGraph, independent_variable::Num)
    sys = ODESystem(Equation[], independent_variable, name = :model)
    return BondGraph(mg, sys)
end

include("OnePorts.jl")
include("Sources.jl")
include("Junctions.jl")
include("TwoPorts.jl")
include("MultiPorts.jl")
include("DerivativeCausality.jl")
include("TransferFunctions.jl")
include("Reactions.jl")
"""

Remove Algebraic Constraint Equations Generated by ModelingToolkit

"""
function remove_algebraic(BG::AbstractBondGraph, model)
    ## Remove Algebraic Constraint Equation
    eqns = equations(model)
    state_var_nodes = filter_vertices(BG.graph, (g, v) -> has_prop(g, v, :state_var)) |> collect
    state_vars = vcat(map(v -> get_prop(BG.graph, v, :state_var), state_var_nodes)...)
    for eqn_index ∈ eachindex(eqns)
        if (eqns[eqn_index].lhs == 0) isa Bool
            vars = get_variables(eqns[eqn_index])
            terms = filter(x -> isa(x, Term), vars)
            # find the extra variables, not a state variable
            guess_loc = 1
            for term_index ∈ eachindex(terms)
                is_state = false
                for sv ∈ state_vars
                    if isa(sv - terms[term_index] == 0, Bool)
                        is_state = true
                    end
                end
                if !is_state
                    res = ModelingToolkit.solve_for(eqns[eqn_index], terms[term_index])
                    eqns = map(eqn -> substitute(eqn, Dict(terms[term_index] => res)), equations(model))
                    model = ODESystem(eqns[1:end-1], name = nameof(model), observed = [observed(model); terms[term_index] ~ res])
                    model = alias_elimination(model)
                    model = structural_simplify(model)
                end
            end
        end
    end
    return model
end

"""

Traverse the BondGraph Structure to create an ODESystem

"""
function graph_to_model(BG::AbstractBondGraph)
    # Find all One Junction Nodes
    one_junctions = filter_vertices(BG.graph, (g, v) -> get_prop(g, v, :type) ∈ [:J1])
    eqns = Equation[]
    for J1 ∈ one_junctions
        out_nodes = outneighbors(BG.graph, J1)
        in_nodes = inneighbors(BG.graph, J1)
        if !isempty(out_nodes)
            out_sum = sum(x -> BG[x].e, out_nodes)
        else
            out_sum = 0
        end
        if !isempty(in_nodes)
            in_sum = sum(x -> BG[x].e, in_nodes)
        else
            in_sum = 0
        end
        push!(eqns, 0 ~ in_sum - out_sum)
        nodes = [out_nodes; in_nodes]
        for i ∈ 2:length(nodes)
            push!(eqns, BG[nodes[i-1]].f ~ BG[nodes[i]].f)
        end
    end
    # Find all Zero Junction Nodes
    zero_junctions = filter_vertices(BG.graph, (g, v) -> get_prop(g, v, :type) ∈ [:J0])
    for J0 ∈ zero_junctions
        out_nodes = outneighbors(BG.graph, J0)
        in_nodes = inneighbors(BG.graph, J0)
        if !isempty(out_nodes)
            out_sum = sum(x -> BG[x].f, out_nodes)
        else
            out_sum = 0
        end
        if !isempty(in_nodes)
            in_sum = sum(x -> BG[x].f, in_nodes)
        else
            in_sum = 0
        end
        push!(eqns, 0 ~ in_sum - out_sum)
        nodes = [out_nodes; in_nodes]
        for i ∈ 2:length(nodes)
            push!(eqns, BG[nodes[i-1]].e ~ BG[nodes[i]].e)
        end
    end
    @named junc_sys = ODESystem(eqns, BG.model.iv, [], [])
    BG.model = extend(BG.model, junc_sys)

    two_ports = filter_vertices(BG.graph, (g, v) -> get_prop(g, v, :type) ∈ [:Re, :TF, :GY, :MTF, :MGY])
    two_ports_sys = map(v -> get_prop(BG.graph, v, :sys), two_ports)
    for sys ∈ two_ports_sys
        BG.model = compose(BG.model, sys)
    end

    element_verts = filter_vertices(BG.graph, (g, v) -> get_prop(g, v, :type) ∈ [:B, :R, :C, :I, :M, :Ce, :Se, :Sf, :MPC, :MPI, :MPR])
    element_sys = map(v -> get_prop(BG.graph, v, :sys), element_verts)
    model = compose(BG.model, element_sys...)
end

"""

Generate an ODE System from the BondGraph Structure

"""
function generate_model!(BG::AbstractBondGraph)
    BG.model = generate_model(BG)
end

function generate_model(BG::AbstractBondGraph)
    model = graph_to_model(BG)
    model = structural_simplify(model)
    while isa(equations(model)[end].lhs == 0, Bool)
        model = remove_algebraic(BG, model)
    end
    return model
end

function generate_model(BG::BioBondGraph)
    model = graph_to_model(BG)
    model = structural_simplify(model)
    RW = SymbolicUtils.Rewriters
    r1 = @acrule log(~x) + log(~y) => log((~x) * (~y))
    r2 = @rule log(~x) - log(~y) => log((~x) / (~y))
    r3 = @rule (~x) * log(~y) => log((~y)^(~x))
    r4 = @rule exp(log(~x)) => ~x
    r5 = @acrule exp((~x) + (~y)) => exp(~x) * exp(~y)
    rw1 = RW.Fixpoint(RW.Chain([r1, r2, r3, r4, r5]))
    rw2 = RW.Prewalk(RW.Chain([r1, r2, r3, r4, r5]))
    rw3 = RW.Postwalk(RW.Chain([r1, r2, r3, r4, r5]))
    eqns = equations(model)
    for i ∈ eachindex(eqns)
        eqns[i] = eqns[i].lhs ~ eqns[i].rhs |> rw3 |> rw2 |> rw1 |> expand
    end
    defaults = model.defaults
    model = ODESystem(eqns, model.iv, states(model), parameters(model), name = nameof(model), observed = observed(model), defaults = defaults)
    return structural_simplify(model)
end

function generate_model!(BG::BioBondGraph)
    BG.model = generate_model(BG)
end
include("IO.jl")

end # module