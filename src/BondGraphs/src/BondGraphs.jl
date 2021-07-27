module BondGraphs

using Base:Real
using Symbolics:Symbolic
using DifferentialEquations
using ModelingToolkit
using SymbolicUtils
using Symbolics
using LinearAlgebra
using LightGraphs
using MetaGraphs
## TODO: Develop better model parameter handling avoid structs of structs....
## TODO: Better interface with LightGraphs nad Flux
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
        get_parameters!,
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
Structure of a bond graph which consists of the ODAE model, (non)linear elements, multi-ports, and initial Conditions
"""
struct Junction
    type::Symbol
    elements::Dict{Symbol,Bool}
    sys::ODESystem
    parameters::Vector{Any}
end

struct Element
    type::Symbol
    sys::ODESystem
    state_var::Vector{Num}
    causality::Bool
end

mutable struct BondGraph
    model::ODESystem
    elements::Dict{Symbol,Element}
    inputs::Vector{Sym{Real,nothing}}
    junctions::Dict{Symbol,Junction}
    initial_state::Dict{Term{Real,Nothing},Number}
    parameters::Vector{Sym}
    graph::MetaGraph{Int64,Float64}
end

Base.getindex(BG::BondGraph, node::Symbol) = get_prop(BG.graph, BG.graph[node, :name], :sys)

## BondGraph constructor
function BondGraph(independent_variable)
    empty_model = ODESystem(Equation[], independent_variable, [], [], systems = []; name = :model)
    mg = MetaGraph(SimpleGraph())
    set_indexing_prop!(mg, :name)
    return BondGraph(empty_model, Dict([]), [], Dict([]), Dict([]), [], mg)
end

include("OnePorts.jl")
include("Sources.jl")
include("Junctions.jl")
include("TwoPorts.jl")
include("MultiPorts.jl")
include("DerivativeCausality.jl")

## Create ODE System From Bond Graph Construction 
function generate_model!(BG::BondGraph)
    junction_0_vertices = filter_vertices(BG.graph, :type, :J0)
    junction_0_sys = map(v -> get_prop(BG.graph, v, :sys), junction_0_vertices)
    junction_1_vertices = filter_vertices(BG.graph, :type, :J1)
    junction_1_sys = map(v -> get_prop(BG.graph, v, :sys), junction_1_vertices)
    TF_vertices = filter_vertices(BG.graph, :type, :TF)
    TF_sys = map(v -> get_prop(BG.graph, v, :sys), TF_vertices)
    GY_vertices = filter_vertices(BG.graph, :type, :GY)
    GY_sys = map(v -> get_prop(BG.graph, v, :sys), GY_vertices)
    MTF_vertices = filter_vertices(BG.graph, :type, :MTF)
    MTF_sys = map(v -> get_prop(BG.graph, v, :sys), MTF_vertices)
    MGY_vertices = filter_vertices(BG.graph, :type, :MGY)
    MGY_sys = map(v -> get_prop(BG.graph, v, :sys), MGY_vertices)
    junction_sys = [junction_1_sys; junction_0_sys; TF_sys; GY_sys; MTF_sys; MGY_sys]
    for sys ∈ junction_sys
        BG.model = extend(sys, BG.model)
    end
    R_vertices = filter_vertices(BG.graph, :type, :R)
    R_sys = map(v -> get_prop(BG.graph, v, :sys), R_vertices)
    C_vertices = filter_vertices(BG.graph, :type, :C)
    C_sys = map(v -> get_prop(BG.graph, v, :sys), C_vertices)
    I_vertices = filter_vertices(BG.graph, :type, :I)
    I_sys = map(v -> get_prop(BG.graph, v, :sys), I_vertices)
    Se_vertices = filter_vertices(BG.graph, :type, :Se)
    Se_sys = map(v -> get_prop(BG.graph, v, :sys), Se_vertices)
    Sf_vertices = filter_vertices(BG.graph, :type, :Sf)
    Sf_sys = map(v -> get_prop(BG.graph, v, :sys), Sf_vertices)

    BG.model = compose(BG.model, [R_sys; C_sys; I_sys; Se_sys; Sf_sys]...)
    display(BG.model)
    nothing
end
## Get parameters
function get_parameters!(BG::BondGraph)
    BG.parameters = map(x -> x => 0.0, parameters(BG.model)) |> Dict
    BG.parameters = keys(BG.parameters) .=> values(BG.parameters)
end
## Simplify Bond Graph System 
simplify_model!(BG::BondGraph) = BG.model = tearing(structural_simplify(BG.model))
# ## Get Independent Variables of system
# function get_states(BG::BondGraph)
#     build_torn_function(BG.model).syms
# end
# ## Set Initial Conditions for Independent Variables
# function set_conditions!(BG::BondGraph, initial_conditions)
#     BG.initial_state = initial_conditions
# end
# ## generate ODEProblem
# function generate_ODE(BG::BondGraph, ps, tspan)
#     ODAEProblem(BG.model, BG.initial_state, tspan)
# end

end # module