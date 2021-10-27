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
export BondGraph,
        add_Bond!, add_R!, add_C!, add_I!, add_M!,
        add_Se!, add_Sf!,
        add_TF!, add_GY!,
        add_MTF!, add_MGY!,
        add_1J!, add_0J!,
        add_C_multiport!, add_I_multiport!,
        generate_model!, generate_model, 
        state_space,
        get_graph, savebondgraph, loadbondgraph

## Function to create a generic Model
mutable struct BondGraph
    graph::MetaDiGraph
    model::ODESystem
end

"""

Get the ODE System Corresponding to the Specific Element

"""
Base.getindex(BG::BondGraph, node::Symbol) = get_prop(BG.graph, BG.graph[node, :name], :sys)

"""

Get the system corresponding to Node - `node`

"""
Base.getindex(BG::BondGraph, node::Int) = get_prop(BG.graph, node, :sys)


"""

Create an empty BondGraph to be populated during the analysis

"""
function BondGraph(independent_variable::Num)
    mg = MetaDiGraph() 
    set_indexing_prop!(mg, :name)
    sys = ODESystem(Equation[], independent_variable, name = :model)
    return BondGraph(mg, sys)
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
"""

Generate an ODE System from the BondGraph Structure

"""
function generate_model!(BG::BondGraph)
    # Find all One Junction Nodes
    one_junctions = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:J1])
    eqns = Equation[]
    for J1 ∈ one_junctions
        out_nodes = outneighbors(BG.graph.graph, J1)
        in_nodes = inneighbors(BG.graph, J1)
        push!(eqns, 0 ~ sum(x->BG[x].e, in_nodes)-sum(x->BG[x].e, out_nodes))
        nodes = [out_nodes; in_nodes]
        for i ∈ 2:length(nodes)
            push!(eqns, BG[nodes[i-1]].f ~ BG[nodes[i]].f)
        end
    end

    # Find all Zero Junction Nodes
    zero_junctions = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:J0])
    for J0 ∈ zero_junctions
        out_nodes = outneighbors(BG.graph, J0)
        in_nodes = inneighbors(BG.graph, J0)
        push!(eqns, 0 ~ sum(x->BG[x].f, in_nodes)-sum(x->BG[x].f, out_nodes))
        nodes = [out_nodes; in_nodes]
        for i ∈ 2:length(nodes)
            push!(eqns, BG[nodes[i-1]].e ~ BG[nodes[i]].e)
        end
    end

    @named junc_sys = ODESystem(eqns, BG.model.iv, [], [])
    BG.model = extend(junc_sys, BG.model)

    two_ports = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:TF, :GY, :MTF, :MGY])
    two_ports_sys = map(v -> get_prop(BG.graph, v, :sys), two_ports)
    for sys ∈ two_ports_sys
        BG.model = extend(sys, BG.model)
    end
    element_verts = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:B, :R, :C, :I, :M, :Se, :Sf, :MPC, :MPI, :MPR])
    element_sys = map(v -> get_prop(BG.graph, v, :sys), element_verts)
    BG.model = compose(BG.model, element_sys...)
    nothing
end

function generate_model(BG::BondGraph)
    # Find all One Junction Nodes
    one_junctions = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:J1])
    eqns = Equation[]
    for J1 ∈ one_junctions
        out_nodes = outneighbors(BG.graph.graph, J1)
        in_nodes = inneighbors(BG.graph, J1)
        push!(eqns, 0 ~ sum(x->BG[x].e, in_nodes)-sum(x->BG[x].e, out_nodes))
        nodes = [out_nodes; in_nodes]
        for i ∈ 2:length(nodes)
            push!(eqns, BG[nodes[i-1]].f ~ BG[nodes[i]].f)
        end
    end

    # Find all Zero Junction Nodes
    zero_junctions = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:J0])
    for J0 ∈ zero_junctions
        out_nodes = outneighbors(BG.graph, J0)
        in_nodes = inneighbors(BG.graph, J0)
        push!(eqns, 0 ~ sum(x->BG[x].f, in_nodes)-sum(x->BG[x].f, out_nodes))
        nodes = [out_nodes; in_nodes]
        for i ∈ 2:length(nodes)
            push!(eqns, BG[nodes[i-1]].e ~ BG[nodes[i]].e)
        end
    end

    @named junc_sys = ODESystem(eqns, BG.model.iv, [], [])
    BG.model = extend(junc_sys, BG.model)

    two_ports = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:TF, :GY, :MTF, :MGY])
    two_ports_sys = map(v -> get_prop(BG.graph, v, :sys), two_ports)
    for sys ∈ two_ports_sys
        BG.model = extend(sys, BG.model)
    end
    element_verts = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:B, :R, :C, :I, :M, :Se, :Sf, :MPC, :MPI, :MPR])
    element_sys = map(v -> get_prop(BG.graph, v, :sys), element_verts)
    compose(BG.model, element_sys...)
end

include("IO.jl")

end # module