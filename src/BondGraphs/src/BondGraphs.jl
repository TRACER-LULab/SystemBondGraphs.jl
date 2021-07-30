module BondGraphs

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
export BondGraph
export add_Bond!
export add_R!
export add_C!
export add_I!
export add_M!
export add_Se!
export add_Sf!
export add_TF!
export add_GY!
export add_MTF!
export add_1J!
export add_0J!
export add_C_multiport!
export add_I_multiport!
export generate_model!
export simplify_model!
export get_parameters!
export get_states
export set_conditions!
export generate_ODE
export check_lhs
export check_rhs
export get_args
export get_implicit
export get_diff
export resolve_derivative_causality!
export transfer_function
export make_graph
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

## BondGraph constructor
function BondGraph(independent_variable)
    empty_model = ODESystem(Equation[], independent_variable, [], [], systems=[])
    mg = MetaGraph(SimpleGraph())
    set_indexing_prop!(mg, :name)
    return BondGraph(empty_model, Dict([]), [], Dict([]), Dict([]), [], mg)
end

## Add Generic Bond to Model
function add_Bond!(BG::BondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    sys = ODESystem(Equation[], BG.model.iv, [e, f], [], name=name)
    BG.elements[name] = Element(:B, sys, [], false)   
    nothing
end

## Add R_element to Model
function add_R!(BG::BondGraph, name; causality=false)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters R
    eqns = [e ~ R * f ]
    sys  = ODESystem(eqns, BG.model.iv, [e, f], [R], name=name)
    BG.elements[name] = Element(:R, sys, [], causality)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end
function add_R!(BG::BondGraph, Φr, params; causality=false)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [e ~ Φr(e, f, BG.model.iv)]
    sys =  ODESystem(eqns, BG.model.iv, [e, f], params, name=name)
    BG.elements[name] = Element(:R, sys, [], causality)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

## Add C-element to Model
function add_C!(BG::BondGraph, name; causality=false)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    @parameters C
    D = Differential(BG.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ q / C
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, q], [C], name=name)
    BG.elements[name] =  Element(:C, sys, [sys.q], causality)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end
function add_C!(BG::BondGraph, Φc, params, name; causality=false)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    D = Differential(BG.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ Φc(e, q, BG.model.iv) # Integral Causality Form
            # 0.0 ~ Φc(e, q, BG.model.iv)
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, q], [], name=name)
    BG.elements[name] = Element(:C, sys, [sys.q], causality)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

## Add I-element to model
function add_I!(BG::BondGraph, name; causality=false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv)
    @parameters I
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            f ~ p / I
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, p], [I], name=name)
    BG.elements[name] = Element(:I, sys, [sys.p], causality)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end
function add_I!(BG::BondGraph, Φi, params, name; causality=false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv)
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            f ~ Φi(p, f, BG.model.iv) # Integral Causality Form
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, p], [], name=name)
    BG.elements[name] = Element(:I, sys, [sys.p], causality)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

## Add M-element to model
function add_M!(BG::BondGraph, name; causality=false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv) q(BG.model.iv)
    @parameters M
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            D(q) ~ f,
            p ~ M * q
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, p, q], [M], name=name)
    BG.elements[name] = Element(:M, sys, [sys.p], causality)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing 
end
function add_M!(BG::BondGraph, Φm, params, name; causality=false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv) q(BG.model.iv)
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            D(q) ~ f,
            p ~ Φi(p, q, BG.model.iv) # Integral Causality Form
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, p, q], [], name=name)
    BG.elements[name] = Element(:M, sys, [sys.p], causality)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

## Add Effort Source
function add_Se!(BG::BondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters Se(BG.model.iv)
    eqns = [0 ~ e - Se]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Se], name=name)
    BG.elements[name] = Element(:Se, sys, [], false)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

function add_Se!(BG::BondGraph, Se::Number, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ e - Se]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Se], name=name)
    BG.elements[name] = Element(:Se, sys, [], false)
    push!(BG.inputs, parameters(sys))
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

function add_Se!(BG::BondGraph, Se, params::Vector{}, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ e - Se(BG.model.iv, params)]
    sys = ODESystem(eqns, BG.model.iv, [e, f], params, name=name)
    BG.elements[name] = Element(:Se, sys, [], false)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

## Add Flow Source
function add_Sf!(BG::BondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters Sf(BG.model.iv)
    eqns = [0 ~ f - Sf]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Sf], name=name)
    BG.elements[name] = Element(:Se, sys, [], false)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

function add_Sf!(BG::BondGraph, Sf::Number, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ f - Sf]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Sf], name=name)
    BG.elements[name] = Element(:Sf, sys, [], false)
    push!(BG.inputs, parameters(sys))
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

function add_Sf!(BG::BondGraph, Sf, params::Vector{}, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ f - Sf(BG.model.iv, params)]
    sys = ODESystem(eqns, BG.model.iv, [e, f], params, name=name)
    BG.elements[name] = Element(:Sf, sys, [], false)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

## Add 1-junction to model
function add_1J!(BG::BondGraph, elements::Dict{Symbol,Bool}, name)
    elems = collect(keys(elements))
    systems = map(x -> BG.elements[x].sys, elems)
    eqns = [
            0 ~ sum(x -> BG.elements[x].sys.e * (-1).^(!elements[x]), keys(elements)) # Sum of all flows is 0
            ]
    for i ∈ 1:length(elems) - 1
        push!(eqns, BG.elements[elems[i]].sys.f ~ BG.elements[elems[i + 1]].sys.f) # effort equality
    end
    sys = ODESystem(eqns, BG.model.iv)
    BG.junctions[name] = Junction(:J1, elements, sys, [])
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    for j ∈ keys(elements)
        add_edge!(BG.graph, BG.graph[name, :name], BG.graph[j, :name])
    end
    nothing
end

## Add 0-junction to model
function add_0J!(BG::BondGraph, elements::Dict{Symbol,Bool}, name)
    eqns = [
            0 ~ sum(x -> BG.elements[x].sys.f * (-1).^(elements[x]), keys(elements)) # Sum of all flows is 0
            ]

    elems = collect(keys(elements))
    for i ∈ 1:length(elems) - 1
        push!(eqns, BG.elements[elems[i]].sys.e ~ BG.elements[elems[i + 1]].sys.e) # effort equality
    end
    # systems = map(x -> BG.elements[x].sys, elems)
    BG.junctions[name] = Junction(:J0, elements, ODESystem(eqns, BG.model.iv, [], [], name=name), [])
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    for j ∈ keys(elements)
        add_edge!(mg, mg[name, :name], mg[j, :name])
    end
    nothing
end

## Add Transformer to Model
function add_TF!(BG::BondGraph, m, elements::Dict{Symbol,Bool}, name)
    # @parameters m
    elems = collect(keys(elements))
    elem_sys = map(x -> BG.elements[x].sys, elems)
    directions = map(x -> elements[x], elems)
    eqns = [
        0 ~ BG.elements[elems[1]].sys.e * (-1)^directions[1] + m * BG.elements[elems[2]].sys.e * (-1)^directions[2], 
        0 ~ m * BG.elements[elems[1]].sys.f * (-1)^directions[1] + BG.elements[elems[2]].sys.f * (-1)^directions[2]
    ]
    element_sys = map(x -> BG.elements[x].sys, elems)
    sys = ODESystem(eqns, BG.model.iv,  [BG.elements[elems[1]].sys.e, BG.elements[elems[2]].sys.f, BG.elements[elems[1]].sys.f, BG.elements[elems[2]].sys.e], [m], name=name)
    BG.junctions[name] = Junction(:TF, elements, sys, parameters(sys))
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

## Add Gyrator to Model Function add_TF(BG, m, elements::Dict{Symbol, Bool},name)
function add_GY!(BG, r, elements::Dict{Symbol,Bool}, name)
    # @parameters r
    elems = collect(keys(elements))
    directions = map(x -> elements[x], elems)
    eqns = [
        0 ~ BG.elements[elems[1]].sys.e * (-1)^directions[1] + (-1) * r * BG.elements[elems[2]].sys.f * (-1)^directions[2], 
        0 ~ r * BG.elements[elems[1]].sys.f * (-1)^directions[1] + (-1) * BG.elements[elems[2]].sys.e * (-1)^directions[2]
    ]
    element_sys = map(x -> BG.elements[x].sys, elems)
    sys = ODESystem(eqns, BG.model.iv, [BG.elements[elems[1]].sys.e, BG.elements[elems[2]].sys.f, BG.elements[elems[1]].sys.f, BG.elements[elems[2]].sys.e], [r], name=name)
    BG.junctions[name] = Junction(:GY, elements, sys, parameters(sys))
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

## Add MTF - Modulated Transformed Modulus
function add_MTF!(BG::BondGraph, m, states, ps, elements, name)
    # @parameters m
    elems = map(x -> x.first, elements)
    elem_sys = map(x -> BG.elements[x].sys, elems)
    directions = map(x -> x.second, elements)
    eqns = [
        0 ~ BG.elements[elems[1]].sys.e * (-1)^directions[1] + m * BG.elements[elems[2]].sys.e * (-1)^directions[2], 
        0 ~ m * BG.elements[elems[1]].sys.f * (-1)^directions[1] + BG.elements[elems[2]].sys.f * (-1)^directions[2]
    ]
    element_sys = map(x -> BG.elements[x].sys, elems)
    sys = ODESystem(eqns, BG.model.iv,  [BG.elements[elems[1]].sys.e, BG.elements[elems[2]].sys.f, BG.elements[elems[1]].sys.f, BG.elements[elems[2]].sys.e, states...], ps, name=name)
    BG.junctions[name] = Junction(:TF, Dict(elems .=> directions), sys, parameters(sys))
    add_vertex!(BG.graph)
    set_prop!(BG.graph, length(BG.graph.graph.fadjlist), :name, name)
    nothing
end

## Create C- multiport
function add_C_multiport!(BG::BondGraph, elements, parameters, name; ϕi=(e, q, params) -> [], ϕk=(e, q, params) -> [])
    # Do the usual setup
    D = Differential(BG.model.iv)
    # Sort Elements 
    𝐪_1j = filter(x -> x.second == false, elements)
    j = length(𝐪_1j)
    𝐞_jp1n = filter(x -> x.second == true, elements)
    n = length(elements)
    # Repack elements based on (7.20) & (7.21)
    elements = [𝐪_1j;𝐞_jp1n]
    # Create variable vectors 
    @variables 𝐪[1:length(elements)](BG.model.iv)
    # @variables 𝐟[1:length(elements)](BG.model.iv)
    # @variables 𝐞[1:length(elements)](BG.model.iv)
    𝐞 = map(i -> BG.elements[elements[i].first].sys.e, eachindex(elements))
    # Create Derivative Relationships for displacement d/dt(q_i) = f_i
    deriv_eqns = map(i -> D(𝐪[i]) ~ BG.elements[elements[i].first].sys.f, eachindex(elements))
    # Create Relationships for (7.20) e_i = ϕ_i(q_1j, e_jn, p)
    𝐞_1j = ϕi(𝐞[j + 1:n], 𝐪[1:j],  parameters)
    e_eqns = map(i -> BG.elements[elements[i].first].sys.e ~ 𝐞_1j[i], 1:j)
    𝐪_jp1n = ϕk(𝐞[j + 1:n], 𝐪[1:j], parameters)
    q_eqns = map(i -> 𝐪[j + i] ~ 𝐪_jp1n[i], 1:n - j)
    eqns = [deriv_eqns; e_eqns; q_eqns]
    eqns = convert(Vector{Equation}, eqns)
    subsys = map(i -> BG.elements[elements[i].first].sys, eachindex(elements))
    sys = compose(ODESystem(eqns, BG.model.iv, collect(𝐪), [], name=name), subsys)
    BG.elements[name] = Element(:C, sys, collect(𝐪), false)
    nothing
end

## Create I-multiport
function add_I_multiport!(BG::BondGraph, elements, parameters; ϕi=(p, f, params) -> [], ϕk=(p, f, params) -> [], name)
    # Do the usual setup
    D = Differential(BG.model.iv)
    # Sort Elements 
    𝐩_1j = filter(x -> x.second == false, elements)
    j = length(𝐩_1j)
    𝐟_jp1n = filter(x -> x.second == true, elements)
    n = length(elements)
    # Repack elements based on (7.20) & (7.21)
    elements = [𝐩_1j;𝐟_jp1n]
    # Create variable vectors 
    @variables 𝐩[1:length(elements)](BG.model.iv)
    𝐟 = map(i -> BG.elements[elements[i].first].sys.f, eachindex(elements))
    # Create Derivative Relationships for displacement d/dt(q_i) = f_i
    deriv_eqns = map(i -> D(𝐩[i]) ~ BG.elements[elements[i].first].sys.e, eachindex(elements))
    # Create Relationships for (7.20) e_i = ϕ_i(q_1j, e_jn, p)
    𝐟_1j = ϕi(𝐩[1:j], 𝐟[j + 1:n], parameters)
    e_eqns = map(i -> BG.elements[elements[i].first].sys.f ~ 𝐟_1j[i], 1:j)
    𝐩_jn = ϕk(𝐩[1:j], 𝐟[j + 1:n], parameters)
    q_eqns = map(i -> 𝐩[j + i] ~ 𝐩_jn[i], 1:n - j)
    return [deriv_eqns; e_eqns; q_eqns]
end

## Create ODE System From Bond Graph Construction 
function generate_model!(BG::BondGraph)
    junctions = collect(values(BG.junctions))
    junc_eqns = reduce(vcat, map(x -> equations(x.sys), junctions))
    junc_sys = reduce(vcat, map(x -> x.sys, junctions))
    
    elem_sys = map(x -> x.sys, collect(values(BG.elements)))

    BG.model = compose(ODESystem(junc_eqns, BG.model.iv, [], []), [elem_sys; junc_sys])
    nothing
end
## Get parameters
function get_parameters!(BG::BondGraph)
    BG.parameters = map(x -> x => 0.0, parameters(BG.model)) |> Dict
    BG.parameters = keys(BG.parameters) .=> values(BG.parameters)
end
## Simplify Bond Graph System 
function simplify_model!(BG::BondGraph)
    BG.model = tearing(structural_simplify(BG.model))
end
## Get Independent Variables of system
function get_states(BG::BondGraph)
    build_torn_function(BG.model).syms
end
## Set Initial Conditions for Independent Variables
function set_conditions!(BG::BondGraph, initial_conditions)
    BG.initial_state = initial_conditions
end
## generate ODEProblem
function generate_ODE(BG::BondGraph, ps, tspan)
    ODAEProblem(BG.model, BG.initial_state, tspan)
end

## Resolve Implicit Equations and Derivative causality

function get_args(term)
    # display("start")
    arg_sub = SymbolicUtils.arguments(term)
    args = []
    for arg ∈ arg_sub
        if !(arg isa Number) && !(arg isa Sym)
            if length(SymbolicUtils.arguments(arg)) > 1
                push!(args, get_args(arg)...)
            else 
                push!(args, arg)
            end
        end
    end
    # display("end")
    return args
end

function check_lhs(eqn, term)
    # LHS Check
    if eqn.lhs isa Term
        if (eqn.lhs - term == 0) isa Bool
            return eqn.lhs ~ eqn.rhs
        elseif (eqn.lhs + term == 0) isa Bool
            return eqn.lhs ~ -1 * eqn.rhs
        end
    end
    return false
end

function check_rhs(eqn, term)
    for arg ∈ get_variables(eqn.rhs)
        if (((arg - term == 0) isa Bool) || ((arg + term == 0) isa Bool))
            return term ~ Symbolics.solve_for(eqn, term)
        end
    end
    return false 
end

function find_eqn(term, eqns)
    eqn = Equation[]
    for i ∈ eachindex(eqns)
        res_lhs = check_lhs(eqns[i], term)
        if res_lhs isa Bool
            res_rhs = check_rhs(eqns[i], term)
            if !(res_rhs isa Bool)
                sub_dict[res_rhs.lhs] = res_rhs.rhs
                break
            end
        else 
            return res_lhs.lhs => res_lhs.rhs
        end
    end
end

function get_diff(eqns)
    diff_eqns = []
    alg_eqns = []
    for eqn ∈ eqns
        if !(eqn.lhs isa Real)
            # Check for differential equations
            if SymbolicUtils.operation(eqn.lhs)  isa Differential
                push!(diff_eqns, eqn)
            # Check for equations with state variables
            else 
                push!(alg_eqns, 0 ~ eqn.lhs - eqn.rhs)
            end
        else 
            push!(alg_eqns, eqn)
        end
    end
    return diff_eqns, alg_eqns
end

function get_implicit(state_vars, srch_eqns, alg_eqns)
    res_eqns = []
    for eqn ∈ srch_eqns
        sub_dict = Dict([])
        old_len = -1.0
        eqns = copy(alg_eqns)
        search_terms = [get_args(eqn.rhs)[1]]
        while length(keys(sub_dict)) != old_len
            old_len = length(keys(sub_dict))
            for search_term ∈ search_terms
                for i ∈ eachindex(eqns)
                    res_lhs = check_lhs(eqns[i], search_term)
                    if res_lhs isa Bool
                        res_rhs = check_rhs(eqns[i], search_term)
                        if !(res_rhs isa Bool)
                            sub_dict[res_rhs.lhs] = res_rhs.rhs
                            popat!(eqns, i)
                            break
                        end
                    else 
                        sub_dict[res_lhs.lhs] = res_lhs.rhs
                        popat!(eqns, i)
                        break
                    end
                end
            end
            eqn = simplify(substitute(eqn, sub_dict))
            eqns = map(x -> expand_derivatives(expand(simplify(substitute(x, sub_dict), expand=true))), eqns)
            args = get_args(eqn.rhs)
            final_args = []
            for arg ∈ args
                check = 0
                for sv ∈ state_vars
                    if !((arg - sv == 0) isa Bool)
                        check += 1
                    end
                end
                (check == length(state_vars)) ? push!(final_args, arg) : ()
            end
        search_terms = union(search_terms, final_args)
        end
        eqns = map(x -> expand_derivatives(expand(simplify(substitute(x, Dict([eqn.lhs => eqn.rhs])),  expand=true))), eqns)
        eqn = expand(simplify(substitute(eqn, sub_dict)))
        push!(res_eqns, eqn)
    end
    return res_eqns
end

function resolve_derivative_causality!(BG::BondGraph)
    # Find the Elements with Derivative Causality aka element.causality=true
    # constiuitive_equations = []
    eqns = equations(BG.model)
    state_vars = reduce(vcat, map(x -> x.state_var, filter(x -> (x.causality == false), collect(values(BG.elements)))))
    D = Differential(BG.model.iv)
    diff_indexes = []
    implicit_eqns = []
    for element ∈ filter(x -> (x.causality == true), collect(values(BG.elements)))
        if element.type == :I
            constiuitive_equation = element.sys.p ~ element.sys.I * element.sys.f
            old_eqn = element.sys.f ~ element.sys.p / element.sys.I
            index_old_eqn = indexin([old_eqn], eqns)[1]
            diff_eqns, alg_eqns = get_diff([eqns[1:(index_old_eqn - 1)];eqns[(index_old_eqn + 1):end]])
            imp_eqn = get_implicit(state_vars, [constiuitive_equation], alg_eqns)
            e_eqn = element.sys.e ~ expand(expand_derivatives(D(imp_eqn[1].rhs)))
            e_eqn = substitute(e_eqn, Dict(map(j -> diff_eqns[j].lhs => diff_eqns[j].rhs, eachindex(diff_eqns))))
            eqns[index_old_eqn] = old_eqn
            diff_eqn_index = indexin([D(element.sys.p) ~ element.sys.e], eqns)[1]
            push!(diff_indexes, diff_eqn_index)
            push!(implicit_eqns, e_eqn)
        elseif element.type == :C
            constiuitive_equation = element.sys.q ~ element.sys.C * element.sys.e
            old_eqn = element.sys.e ~ element.sys.q / element.sys.C
            index_old_eqn = indexin([old_eqn], eqns)[1]
            diff_eqns, alg_eqns = get_diff([eqns[1:(index_old_eqn - 1)];eqns[(index_old_eqn + 1):end]])
            imp_eqn = get_implicit(state_vars, [constiuitive_equation], alg_eqns)
            f_eqn = element.sys.f ~ expand(expand_derivatives(D(imp_eqn[1].rhs)))
            f_eqn = substitute(f_eqn, Dict(map(j -> diff_eqns[j].lhs => diff_eqns[j].rhs, eachindex(diff_eqns))))
            eqns[index_old_eqn] = old_eqn
            diff_eqn_index = indexin([D(element.sys.q) ~ element.sys.f], eqns)[1]
            push!(diff_indexes, diff_eqn_index)
            push!(implicit_eqns, f_eqn)
        end
    end
    map(i -> eqns[diff_indexes[i]] = implicit_eqns[i], eachindex(implicit_eqns))
    BG.model = ODESystem(eqns, BG.model.iv, states(BG.model), parameters(BG.model))
end

## Start Implementing LTI - Analysis on form ẋ = Ax+Bu, y = Cx +Du
function transfer_function(BG::BondGraph, ps, C, D)
    sts = states(BG.model)
    eqns = equations(BG.model)
    ins = BG.inputs
    @show ins
    in_dict = Dict(ins .=> 0.0)
    st_dict = Dict(sts .=> 0.0)
    A = []
    B = []
    eqns = map(eqn -> substitute(eqn, ps), eqns)
    for i ∈ eachindex(eqns)
        Aᵢ = []
        for j ∈ eachindex(sts)
            st_dict[sts[j]] = 1.0
            push!(Aᵢ, substitute(substitute(eqns[i].rhs, st_dict), in_dict) |> simplify)
            st_dict[sts[j]] = 0.0
        end
        Bᵢ = []
        for j ∈ eachindex(ins)
            in_dict[ins[j]] = 1.0
            push!(Bᵢ, substitute(substitute(eqns[i].rhs, st_dict), in_dict) |> simplify)
            in_dict[ins[j]] = 0.0
        end
        push!(A, Aᵢ) 
    push!(B, Bᵢ)
    end

    A = reduce(hcat, A) |> permutedims
    B = reduce(hcat, B) |> permutedims
    @variables s
    build_function(C * (s * I(length(sts)) - A)^(-1) * B + D, s)
end
## Graph Theory Implementations
end # module