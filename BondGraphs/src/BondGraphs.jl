module BondGraphs

using DifferentialEquations
using ModelingToolkit
using SymbolicUtils
using Symbolics

export BondGraph
export add_Bond!
export add_R!
export add_C!
export add_I!
# export add_M!
export add_Se!
export add_Sf!
export add_TF!
export add_GY!
export add_1J!
export add_0J!
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
## Function to create a generic Model
"""
Structure of a bond graph which consists of the ODAE model, (non)linear elements, multi-ports, and initial Conditions
"""
struct Junction
    elements::Vector{Symbol}
    directions::Vector{Bool}
    eqns::Vector{Equation}
    parameters::Vector{Any} 
end

struct Element
    type::Symbol
    name::Symbol
    sys::ODESystem
    parameters
    state_var::Vector{Num}
end

mutable struct BondGraph
    model::ODESystem
    elements::Dict{Symbol,Element}
    junctions::Dict{Symbol,Junction}
    initial_state::Dict{Term{Real,Nothing},Number}
    parameters::Vector{Sym}
end


# """
# Creation of an empty bond graph created with the specification of an independent variable

# # Arguments
# - `independent_variable::Symbol`: The independent variable for which the system is solved w.r.t
# """
function BondGraph(independent_variable)
    empty_model = ODESystem(Equation[], independent_variable, [], [], systems=[])
    return BondGraph(empty_model, Dict([]), Dict([]), Dict([]), [])
end

## Add Generic Bond to Model
function add_Bond!(BondGraph, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    sys = ODESystem(Equation[], BondGraph.model.iv, [e, f], [], name=name)
    BondGraph.elements[name] = Element(:B, name, sys, [], [])
end

## Add R_element to Model
"""
add_R! with two arguements assumes that the element is linear
"""
# # function add_R!(BondGraph, R::Number, name)
# #     val = R
# #     @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
# #     @parameters R
# #     eqns = [e ~ R * f]
# #     sys = ODESystem(eqns, BondGraph.model.iv, [e, f], [R], defaults=Dict(R => val), name=name)
# #     # element = Element(name, sys, [R], [])
# #     BondGraph.elements[name] = Element(:R, name, sys, [R], [])
# # end
function add_R!(BondGraph, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    @parameters R
    eqns = [e ~ R * f ]
    sys  = ODESystem(eqns, BondGraph.model.iv, [e, f], [R], name=name)
    BondGraph.elements[name] = Element(:R, name, sys, [R], [])
end
function add_R!(BondGraph, Φr, params, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [e ~ Φr(e, f, BondGraph.model.iv)]
    sys =  ODESystem(eqns, BondGraph.model.iv, [e, f], [], name=name)
    BondGraph.elements[name] = Element(:R, name, sys, params, [])
end

# ## Add C-element to Model
# # function add_C!(BondGraph, C::Number, name)
# #     val = C
# #     @variables e(BondGraph.model.iv) f(BondGraph.model.iv) q(BondGraph.model.iv)
# #     @parameters C
# #     D = Differential(BondGraph.model.iv)
# #     eqns = [
# #             D(q) ~ f,
# #             e ~ q / C
# #             ]
# #     sys = ODESystem(eqns, BondGraph.model.iv, [e, f, q], [C], defaults=Dict(C => val), name=name)
# #     BondGraph.elements[name] = Element(:C, name, sys, [C], [sys.q])
# #     push!(BondGraph.state_vars, BondGraph.elements[name].q)
# # end
function add_C!(BondGraph, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) q(BondGraph.model.iv)
    @parameters C
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ q / C
            ]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f, q], [C], name=name)
    BondGraph.elements[name] =  Element(:C, name, sys, [C], [sys.q])
    # push!(BondGraph.state_vars, BondGraph.elements[name].q)
end
function add_C!(BondGraph, Φc, params, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) q(BondGraph.model.iv)
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ Φc(e, q, BondGraph.model.iv) # Integral Causality Form
            # 0.0 ~ Φc(e, q, BondGraph.model.iv)
            ]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f, q], [], name=name)
    BondGraph.elements[name] = Element(:C, name, sys, params, [sys.q])
    # push!(BondGraph.state_vars, BondGraph.elements[name].q)
end

# ## Add I-element to model
# # function add_I!(BondGraph, I::Number, name)
# #     val = I
# #     @variables e(BondGraph.model.iv) f(BondGraph.model.iv) p(BondGraph.model.iv)
# #     @parameters I
# #     D = Differential(BondGraph.model.iv)
# #     eqns = [
# #             D(p) ~ e,
# #             f ~ p / I
# #             ]
# #     BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f, p], [I], defaults=Dict(I => val), name=name)
# #     push!(BondGraph.state_vars, BondGraph.elements[name].p)
# # end
function add_I!(BondGraph, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) p(BondGraph.model.iv)
    @parameters I
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(p) ~ e,
            f ~ p / I
            ]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f, p], [I], name=name)
    BondGraph.elements[name] = Element(:I, name, sys, [I], [sys.p])
    # push!(BondGraph.state_vars, BondGraph.elements[name].p)
end
function add_I!(BondGraph, Φi, params, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) p(BondGraph.model.iv)
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(p) ~ e,
            f ~ Φi(p, f, BondGraph.model.iv) # Integral Causality Form
            ]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f, p], [], name=name)
    BondGraph.elements[name] = Element(:I, name, sys, params, [sys.p])
    # push!(BondGraph.state_vars, BondGraph.elements[name].p)
end

## Add Memrsistive element
function add_M!(BondGraph, M::Number, name)
    val = M
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) p(BondGraph.model.iv) q(BondGraph.model.iv)
    @parameters M
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(p) ~ e,
            D(q) ~ f,
            0.0 ~ p - M * q
            ]
    BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f, p, q], [M], name=name)
    push!(BondGraph.state_vars, BondGraph.elements[name].p)
    push!(BondGraph.state_vars, BondGraph.elements[name].q)
end
function add_M!(BondGraph, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) p(BondGraph.model.iv) q(BondGraph.model.iv)
    @parameters M
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(p) ~ e,
            D(q) ~ f,
            0.0 ~ p - M * q
            ]
    BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f, p, q], [M], name=name)
    push!(BondGraph.state_vars, BondGraph.elements[name].p)
    push!(BondGraph.state_vars, BondGraph.elements[name].q)
end
function add_M!(BondGraph, Φm, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) p(BondGraph.model.iv) q(BondGraph.model.iv)
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(p) ~ e,
            D(q) ~ f,
            0.0 ~ p - Φm(p, q, BondGraph.model.iv)
            ]
    BondGraph.elements[name] = ODESystem(eqns, BondGraph.model.iv, [e, f, p, q], [], name=name)
    push!(BondGraph.state_vars, BondGraph.elements[name].p)
    push!(BondGraph.state_vars, BondGraph.elements[name].q)
end

## Add Effort Source
function add_Se!(BondGraph, Se::Number, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [0.0 ~ e - Se]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], [], name=name)
    BondGraph.elements[name] = Element(:Se, name, sys, [], [])
end
function add_Se!(BondGraph, Se, params, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [0.0 ~ e - Se(e, f, BondGraph.model.iv, params)]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], [], name=name)
    BondGraph.elements[name] = Element(:Se, name, sys, params, [])
end

## Add Flow Source
function add_Sf!(BondGraph, Sf::Number, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [0.0 ~ f - Sf]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], [], name=name)
    BondGraph.elements[name] = Element(:Sf, name, sys, [], [])
end
function add_Sf!(BondGraph, Sf, ps::Vector, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [0.0 ~ f - Sf(e, f, BondGraph.model.iv)]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], [], name=name)
    BondGraph.elements[name] = Element(:Sf, name, sys, ps, [])
end

## Add 1-junction to model
function add_1J!(BondGraph, elements::Dict{Symbol,Bool}, name)
    eqns = [
            0.0 ~ sum(x -> BondGraph.elements[x].sys.e * (-1).^(elements[x]), keys(elements)) # Sum of all flows is 0
            ]
    elems = collect(keys(elements))
    for i ∈ 1:length(elems) - 1
        push!(eqns, BondGraph.elements[elems[i]].sys.f ~ BondGraph.elements[elems[i + 1]].sys.f) # effort equality
    end
    BondGraph.junctions[name] = Junction(collect(keys(elements)), collect(values(elements)), eqns, [])
end

## Add 0-junction to model
function add_0J!(BondGraph, elements::Dict{Symbol,Bool}, name)
    eqns = [
            0.0 ~ sum(x -> BondGraph.elements[x].sys.f * (-1).^(elements[x]), keys(elements)) # Sum of all flows is 0
            ]

    elems = collect(keys(elements))
    for i ∈ 1:length(elems) - 1
        push!(eqns, BondGraph.elements[elems[i]].sys.e ~ BondGraph.elements[elems[i + 1]].sys.e) # effort equality
    end
    BondGraph.junctions[name] = Junction(collect(keys(elements)), collect(values(elements)), eqns, [])
end

## Add Transformer to Model
function add_TF!(BondGraph, m, elements::Dict{Symbol,Bool}, name)
    elems = collect(keys(elements))
    directions = map(x -> elements[x], elems)
    eqns = [
        0.0 ~ BondGraph.elements[elems[1]].sys.e * (-1)^directions[1] - m * BondGraph.elements[elems[2]].sys.e * (-1)^directions[2], 
        0.0 ~ m * BondGraph.elements[elems[1]].sys.f * (-1)^directions[1] - BondGraph.elements[elems[2]].sys.f * (-1)^directions[2]
    ]
    BondGraph.junctions[name] = Junction(collect(keys(elements)), collect(values(elements)), eqns, [m])
end

# Add Gyrator to Model Function add_TF(BondGraph, m, elements::Dict{Symbol, Bool},name)
function add_GY!(BondGraph, r, elements::Dict{Symbol,Bool}, name)
    elems = collect(keys(elements))
    directions = map(x -> elements[x], elems)
    eqns = [
        0.0 ~ BondGraph.elements[elems[1]].sys.e * (-1)^directions[1] + (-1) * r * BondGraph.elements[elems[2]].sys.f * (-1)^directions[2], 
        0.0 ~ r * BondGraph.elements[elems[1]].sys.f * (-1)^directions[1] + (-1) * BondGraph.elements[elems[2]].sys.e * (-1)^directions[2]
    ]
    BondGraph.junctions[name] = Junction(collect(keys(elements)), collect(values(elements)), eqns, [r])
end

## Create ODE System From Bond Graph Construction 
function generate_model!(BondGraph)
    elements = collect(values(BondGraph.elements))
    elements = map(x -> x.sys, elements)
    junctions = collect(values(BondGraph.junctions))
    junctions = reduce(vcat, map(x -> x.eqns, junctions))
    BondGraph.model = ODESystem(junctions, BondGraph.model.iv, systems=elements)
end
## Get parameters
function get_parameters!(BondGraph)
    BondGraph.parameters = map(x -> x => 0.0, parameters(BondGraph.model)) |> Dict
    BondGraph.parameters = keys(BondGraph.parameters) .=> values(BondGraph.parameters)
end
## Simplify Bond Graph System 
function simplify_model!(BondGraph)
    BondGraph.model = structural_simplify(BondGraph.model)
end
## Get Independent Variables of system
function get_states(BondGraph)
    build_torn_function(BondGraph.model).syms
end
## Set Initial Conditions for Independent Variables
function set_conditions!(BondGraph, initial_conditions)
    BondGraph.initial_state = initial_conditions
end
## generate ODEProblem
function generate_ODE(BondGraph, ps, tspan)
    ODAEProblem(BondGraph.model, BondGraph.initial_state, tspan)
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
    for arg ∈ SymbolicUtils.arguments(eqn.rhs)
        if ((arg - term == 0) isa Bool) || ((arg + term == 0) isa     Bool)
            return term ~ -Symbolics.solve_for(eqn, term)
        end
    end
    return false 
end

function get_args(term)
    arg_sub = SymbolicUtils.arguments(term)
    args = []
    for arg ∈ arg_sub

        if !(arg isa Number) && !(arg isa Sym)
            if length(SymbolicUtils.arguments(arg)) > 1
                get_args(arg)
            else 
                push!(args, arg)
            end
        end
    end
    return args
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

function get_implicit(state_vars, dc_eqn, eqns)
    sub_dict = Dict([])
    old_len = -1.0
    search_terms = [get_args(dc_eqn.rhs)[1]]
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
        dc_eqn = simplify(substitute(dc_eqn, sub_dict), expand=true)
        eqns = map(x -> expand_derivatives(expand(simplify(substitute(x, sub_dict), expand=true))), eqns)
        args = get_args(dc_eqn.rhs)
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
    eqns = map(x -> expand_derivatives(expand(simplify(substitute(x, Dict([dc_eqn.lhs => dc_eqn.rhs])), expand=true))), eqns)
    dc_eqn = simplify(substitute(dc_eqn, sub_dict), expand=true)
    return dc_eqn, eqns, sub_dict
end

end # module