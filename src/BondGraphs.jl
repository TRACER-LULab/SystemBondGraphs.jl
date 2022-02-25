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
    generate_model!, generate_model,
    state_space,
    get_graph, savebondgraph, loadbondgraph,
    remove_algebraic, remove_casuality,
    derivative_casuality

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

include("LinearOnePorts.jl")
include("NonlinearOnePorts.jl")
include("Sources.jl")
include("Junctions.jl")
include("TwoPorts.jl")
include("MultiPorts.jl")
include("TransferFunctions.jl")
include("Reactions.jl")
include("ModelGeneration.jl")
include("IO.jl")
include("DerivativeCausality.jl")
include("NewDC.jl")
end # module