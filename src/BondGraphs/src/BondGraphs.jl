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
mutable struct BondGraph
    graph::MetaDiGraph
    model::ODESystem
end

"""

Get the ODE System Corresponding to the Specific Element

"""
Base.getindex(BG::BondGraph, node::Symbol) = get_prop(BG.graph, BG.graph[node, :name], :sys)

"""

Create an empty BondGraph to be populated during the analysis

"""
function BondGraph(independent_variable)
    mg = MetaDiGraph()
    set_indexing_prop!(mg, :name)
    sys = ODESystem(Equation[], independent_variable, name = :model)
    return BondGraph(mg, sys)
end

function BondGraph(mg::AbstractMetaGraph, independent_variable)
    sys = ODESystem(Equation[], independent_variable, name = :model)
    nlabels = [string(mg.vprops[i][:name]) for i ∈ 1:length(keys(mg.vprops))]
    junctions = filter(x->(x[1:2] == "J0" || x[1:2] == "J1"), nlabels)
    one_junctions = filter(x->x[2] == "1", junctions)
    zero_junctions = filter(x->x[2] == "0", junctions)
    

    return BondGraph(mg, sys)
end

include("OnePorts.jl")
include("Sources.jl")
include("Junctions.jl")
include("TwoPorts.jl")
include("MultiPorts.jl")
include("DerivativeCausality.jl")

"""

Generate an ODE System from the BondGraph Structure

"""
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

include("Import.jl")

end # module