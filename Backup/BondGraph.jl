using DifferentialEquations
using Plots
using RecursiveArrayTools
using Symbolics
using SymbolicUtils
# using GLMakie
using ModelingToolkit
## Function to create a generic Model
mutable struct BondGraph
    model::ODESystem
    elements::Dict{Symbol, ODESystem}
    junctions::Dict{Symbol, Vector{Equation}}
    initial_state::Dict{Term{Real, Nothing}, Number}
end

function BondGraph(independent_variable)
    empty_model = ODESystem(Equation[], independent_variable, [], [], systems = [])
    return BondGraph(empty_model,Dict([]), Dict([]), Dict([]))
end

## Add Generic Bond to Model
function add_Bond!(BondGraph, name)
    @variables e(t) f(t)
    BondGraph.elements[name] = ODESystem(Equation[], t, [e, f], [], defaults=[e => 0.0, f => 0.0], name=name)
end 

## Add R_element to Model
function add_R!(BondGraph, R::Number, name)
    val = R
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    @parameters R
    eqns = [0.0 ~ R * f - e]
    BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f], [R], defaults=Dict(R => val, e => 0.0, f => 0.0), name = name)
end
function add_R!(BondGraph, Φr, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [0.0 ~ Φr(f, t) - e]
    BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f], [], defaults=Dict(e => 0.0, f => 0.0), name = name)
end

## Add C-element to Model
function add_C!(BondGraph, C::Number, name)
    val = C
    @variables e(t) f(t) q(t)
    @parameters C
    D = Differential(t)
    eqns = [
            D(q) ~ f,
            0.0 ~ q-C*e
            ]
    BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f, q], [C], defaults=Dict(C => val, e => 0.0, f => 0.0, q=>0.0), name = name)
end
function add_C!(BondGraph, Φc, name)
    @variables e(t) f(t) q(t)
    D = Differential(t)
    eqns = [
            D(q) ~ f,
            0.0 ~ q-Φc(q, e, t)
            ]
    BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f, q], [], defaults=Dict(e => 0.0, f => 0.0, q=>0.0), name = name)
end

## Add I-element to model
function add_I!(BondGraph, I::Number, name)
    val = I
    @variables e(t) f(t) p(t)
    @parameters I
    D = Differential(t)
    eqns = [
            D(p) ~ e,
            0.0 ~ p-f*I
            ]
    BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f, p], [I], defaults=Dict(I => val, e => 0.0, f => 0.0, p => 0.0), name=name)
end

function add_I!(BondGraph, Φi, name)
    @variables e(t) f(t) p(t)
    D = Differential(t)
    eqns = [
            D(p) ~ e,
            0.0 ~ p-Φi(p, f, t)
            ]
    BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f, p], [], defaults=Dict(e => 0.0, f => 0.0, p => 0.0), name=name)
end

## Add Effort Source
function add_Se!(BondGraph, Se, name)
    @variables e(t) f(t)
    eqns = [0.0 ~ e - Se]
    BondGraph.elements[name] = ODESystem(eqns, t, [e, f], [], name=name)
end

## Add Flow Source
function add_Sf!(BondGraph, Sf, name)
    @variables e(t) f(t)
    eqns = [0.0 ~ f - Sf]
    BondGraph.elements[name] = ODESystem(eqns, t, [e, f], [], name=name)
end

## Add 1-junction to model
function add_1J!(BondGraph, elements::Dict{Symbol, Bool}, name)
    eqns = [
            0.0 ~ sum(x -> BondGraph.elements[x].e * (-1).^(elements[x]), keys(elements)) # Sum of all flows is 0
            ]

    elems = collect(keys(elements))
    for i ∈ 1:length(elems) - 1
        push!(eqns, 0.0 ~ BondGraph.elements[elems[i]].f - BondGraph.elements[elems[i + 1]].f) # effort equality
    end
    BondGraph.junctions[name] = eqns
end

## Add 0-junction to model
function add_0J!(BondGraph, elements::Dict{Symbol, Bool}, name)
    eqns = [
            0.0 ~ sum(x -> BondGraph.elements[x].f * (-1).^(elements[x]), keys(elements)) # Sum of all flows is 0
            ]

    elems = collect(keys(elements))
    for i ∈ 1:length(elems) - 1
        push!(eqns, 0.0 ~ BondGraph.elements[elems[i]].e - BondGraph.elements[elems[i + 1]].e) # effort equality
    end
    BondGraph.junctions[name] = eqns
end

## Add Transformer to Model
function add_TF!(BondGraph, m, elements::Dict{Symbol, Bool},name)
    elems = collect(keys(elements))
    directions = map(x->elements[x], elems)
    eqns = [
        0 ~ BondGraph.elements[elems[1]].e * (-1)^directions[1] - m * BondGraph.elements[elems[2]].e * (-1)^directions[2], 
        0 ~ m * BondGraph.elements[elems[1]].f * (-1)^directions[1] - BondGraph.elements[elems[2]].f * (-1)^directions[2]
    ]
    BondGraph.junctions[name] = eqns
end

## Add Gyrator to Model Function add_TF(BondGraph, m, elements::Dict{Symbol, Bool},name)
function add_GY!(BondGraph, r, elements::Dict{Symbol, Bool},name)
    elems = collect(keys(elements))
    directions = map(x->elements[x], elems)
    eqns = [
        0 ~ BondGraph.elements[elems[1]].e * (-1)^directions[1] - r * BondGraph.elements[elems[2]].f * (-1)^directions[2], 
        0 ~ r * BondGraph.elements[elems[1]].f * (-1)^directions[1] - BondGraph.elements[elems[2]].e * (-1)^directions[2]
    ]
    BondGraph.junctions[name] = eqns
end

## Create ODE System From Bond Graph Construction 
function generate_model!(BondGraph)
    elements = collect(values(BondGraph.elements))
    junctions = collect(values(BondGraph.junctions))
    junctions = reduce(vcat, junctions)
    BondGraph.model = ODESystem(junctions, BondGraph.model.iv, systems = elements)
end
function generate_model!(BondGraph, p)
    elements = collect(values(BondGraph.elements))
    junctions = collect(values(BondGraph.junctions))
    junctions = reduce(vcat, junctions)
    BondGraph.model = ODESystem(junctions, BondGraph.model.iv, ps = p, systems = elements)
end
## Simplify Bond Graph System 
function simplify_model!(BondGraph)
    BondGraph.model = structural_simplify(BondGraph.model)
end
## Get Independent Variables of system
function get_states(BondGraph)
    states(BondGraph.model)
end
## Set Initial Conditions for Independent Variables
function set_conditions!(BondGraph, initial_conditions)
    BondGraph.initial_state = initial_conditions
end
## generate ODEProblem
function generate_ODE(BondGraph, tspan)
    ODAEProblem(BondGraph.model, BondGraph.initial_state, tspan)
end