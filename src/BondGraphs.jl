module BondGraphs

using DifferentialEquations
using ModelingToolkit
using SymbolicUtils
using Symbolics
using LinearAlgebra
using Graphs
using MetaGraphs
using FileIO

export BondGraph, BioBondGraph
export add_Bond!, add_R!, add_C!, add_I!, add_M!
export add_Se!, add_Sf!
export add_TF!, add_GY!
export add_MTF!, add_MGY!
export add_1J!, add_0J!
export add_IP!
export generate_model!, generate_model
export state_space
export get_graph, savebondgraph, loadbondgraph
export remove_algebraic, remove_casuality
export derivative_casuality

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
function BondGraph(independent_variable::Num, name)
    mg = MetaDiGraph()
    set_indexing_prop!(mg, :name)
    sys = ODESystem(Equation[], independent_variable, name=Symbol(name))
    return BondGraph(mg, sys)
end

function BioBondGraph(independent_variable::Num, name::Symbol; R=1.0, T=1.0)
    mg = MetaDiGraph()
    set_indexing_prop!(mg, :name)
    sys = ODESystem(Equation[], independent_variable, name=Symbol(name))
    return BioBondGraph(mg, sys, R, T)
end

"""

Create a BondGraph provided a directed metagraph with node `:name` and `:type` defined for each node

"""
function BondGraph(mg::AbstractMetaGraph, independent_variable::Num, name)
    sys = ODESystem(Equation[], independent_variable, name=Symbol(name))
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
include("InformationProcessor.jl")
end # module
