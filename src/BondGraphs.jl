module BondGraphs

using Reexport
@reexport using ModelingToolkit
# @reexport using Symbolics
@reexport using Graphs
# using SymbolicUtils
using LinearAlgebra
@reexport using MetaGraphsNext
# using FileIO
# using ProtoStructs

export BondGraph, ChemBondGraph
export add_R!, add_C!, add_I!, add_M!
export add_Se!, add_Sf!
export add_TF!, add_GY!
export add_MTF!, add_MGY!
export add_1J!, add_0J!
export add_bond!
export generate_model
# export add_IP!
# export generate_model!
# export state_space
# export get_graph, savebondgraph, loadbondgraph
# export remove_algebraic, remove_casuality
# export derivative_casuality
# export graph_to_model_new
## Function to create a generic Model

struct BondGraphData
    iv::Number
end

struct ChemBondGraphData
    iv::Number
    R::Number
    T::Number
end

struct BondGraphNode
    model::ODESystem
    type::Symbol
end

struct BondGraphEdge #<: Graphs.AbstractEdge
    name::Symbol
    model::ODESystem
end
"""

Get the ODE System Corresponding to the Specific Element

"""
# Base.getindex(bg::AbstractBondGraph, node::Symbol) = bg.graph[node].model
# Base.getindex(bg::AbstractBondGraph, node1::Symbol, node2::Symbol) = bg.graph[node1, node2].model
# Base.copy(bg::BondGraph) = BondGraph(copy(bg.graph), deepcopy(bg.model))
"""

Create an empty BondGraph to be populated during the analysis

"""
function BondGraph(iv::Num)
    mg = MetaGraph(
        DiGraph(),
        label_type=Symbol,
        vertex_data_type=BondGraphNode,
        edge_data_type=BondGraphEdge,
        graph_data=BondGraphData(iv))
    # set_indexing_prop!(mg, :name)
    # sys = ODESystem(Equation[], independent_variable, name=Symbol(name))
    return mg #BondGraph(mg, sys)
end

function ChemBondGraph(iv::Num; R=1.0, T=1.0)
    mg = MetaGraph(
        DiGraph(),
        label_type=Symbol,
        vertex_data_type=BondGraphNode,
        edge_data_type=BondGraphEdge,
        graph_data=ChemBondGraphData(iv, R, T))
    # set_indexing_prop!(mg, :name)
    # sys = ODESystem(Equation[], independent_variable, name=Symbol(name))
    return mg# BioBondGraph(mg, sys, R, T)
end

"""

Create a BondGraph provided a directed metagraph with node `:name` and `:type` defined for each node

"""
# function BondGraph(mg::AbstractGraph, independent_variable::Num, name)
#     sys = ODESystem(Equation[], independent_variable, name=Symbol(name))
#     return BondGraph(mg, sys)
# end

include("LinearOnePorts.jl")
include("NonlinearOnePorts.jl")
include("Sources.jl")
include("Junctions.jl")
include("TwoPorts.jl")
# include("MultiPorts.jl")
# include("TransferFunctions.jl")
include("Reactions.jl")
include("ModelGeneration.jl")
# include("IO.jl")
# include("DerivativeCausality.jl")
# include("NewDC.jl")
# include("InformationProcessor.jl")
include("Connect.jl")
end # module
