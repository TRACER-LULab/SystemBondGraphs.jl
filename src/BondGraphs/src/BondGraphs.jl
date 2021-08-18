module BondGraphs

using DifferentialEquations
using ModelingToolkit
using SymbolicUtils
using Symbolics
using LinearAlgebra
using LightGraphs
using MetaGraphs
using FileIO
## TODO: Develop better model parameter handling avoid structs of structs....
## TODO: Better interface with LightGraphs and Flux
export BondGraph,
        add_Bond!,
        add_R!,
        add_C!,
        add_I!,
        add_M!,
        add_Se!,
        add_Sf!,
        add_TF!,
        add_GY!,
        add_MTF!,
        add_1J!,
        add_0J!,
        add_C_multiport!,
        add_I_multiport!,
        generate_model!,
        simplify_model!,
        # get_parameters!,
        get_states,
        set_conditions!,
        generate_ODE,
        check_lhs,
        check_rhs,
        get_args,
        get_implicit,
        get_diff,
        resolve_derivative_causality!,
        transfer_function
## Function to create a generic Model
"""
Structure of a bond graph which consists of the system model and a graph representation of the system
"""
mutable struct BondGraph
    graph::MetaGraph
    model::ODESystem
end

Base.getindex(BG::BondGraph, node::Symbol) = get_prop(BG.graph, BG.graph[node, :name], :sys)

## BondGraph constructor
function BondGraph(independent_variable)
    mg = MetaGraph()
    set_indexing_prop!(mg, :name)
    sys = ODESystem(Equation[], independent_variable, name = :model)
    return BondGraph(mg, sys)
end

include("OnePorts.jl")
include("Sources.jl")
include("Junctions.jl")
include("TwoPorts.jl")
include("MultiPorts.jl")
include("DerivativeCausality.jl")

## Create ODE System From Bond Graph Construction 
function generate_model!(BG::BondGraph)
    junction_verts = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:J0, :J1, :TF, :GY, :MTF, :MGY])
    junction_sys = map(v -> get_prop(BG.graph, v, :sys), junction_verts)
    for sys ∈ junction_sys
        BG.model = extend(sys, BG.model)
    end
    element_verts = filter_vertices(BG.graph,  (g, v) -> get_prop(g, v, :type) ∈ [:B, :R, :C, :I, :M, :Se, :Sf, :MPC, :MPI, :MPR])
    element_sys = map(v -> get_prop(BG.graph, v, :sys), element_verts)
    BG.model = compose(BG.model, element_sys...)
    nothing
end

## Simplify Bond Graph System 
simplify_model!(BG::BondGraph) = BG.model = tearing(structural_simplify(BG.model))

include("Import.jl")

end # module