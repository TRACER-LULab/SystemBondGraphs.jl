module BondGraphs

using Symbolics:Symbolic
using DifferentialEquations
using ModelingToolkit
using SymbolicUtils
using Symbolics
using LinearAlgebra
## TODO: Develop better model parameter handling avoid structs of structs....
## TODO: Better interface with LightGraphs nad Flux
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
export add_MTF!
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
export get_diff
export resolve_derivative_causality!
export transfer_function
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
end


# """
# Creation of an empty bond graph created with the specification of an independent variable

# # Arguments
# - `independent_variable::Symbol`: The independent variable for which the system is solved w.r.t
# """
function BondGraph(independent_variable)
    empty_model = ODESystem(Equation[], independent_variable, [], [], systems=[])
    return BondGraph(empty_model, Dict([]), [], Dict([]), Dict([]), [])
end

## Add Generic Bond to Model
function add_Bond!(BondGraph, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    sys = ODESystem(Equation[], BondGraph.model.iv, [e, f], [], name=name)
    BondGraph.elements[name] = Element(:B, sys, [], false)
end

## Add R_element to Model
"""
add_R! with two arguments assumes that the element is linear
"""
function add_R!(BondGraph, name; causality=false)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    @parameters R
    eqns = [e ~ R * f ]
    sys  = ODESystem(eqns, BondGraph.model.iv, [e, f], [R], name=name)
    BondGraph.elements[name] = Element(:R, sys, [], causality)
end
function add_R!(BondGraph, Φr, params; causality=false)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [e ~ Φr(e, f, BondGraph.model.iv)]
    sys =  ODESystem(eqns, BondGraph.model.iv, [e, f], params, name=name)
    BondGraph.elements[name] = Element(:R, sys, [], causality)
end

# ## Add C-element to Model
function add_C!(BondGraph, name; causality=false)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) q(BondGraph.model.iv)
    @parameters C
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ q / C
            ]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f, q], [C], name=name)
    BondGraph.elements[name] =  Element(:C, sys, [sys.q], causalitys)
    # push!(BondGraph.state_vars, BondGraph.elements[name].q)
end
function add_C!(BondGraph, Φc, params, name; causality=false)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) q(BondGraph.model.iv)
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ Φc(e, q, BondGraph.model.iv) # Integral Causality Form
            # 0.0 ~ Φc(e, q, BondGraph.model.iv)
            ]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f, q], [], name=name)
    BondGraph.elements[name] = Element(:C, sys, [sys.q], causality)
    # push!(BondGraph.state_vars, BondGraph.elements[name].q)
end

# ## Add I-element to model
function add_I!(BondGraph, name; causality=false)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) p(BondGraph.model.iv)
    @parameters I
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(p) ~ e,
            f ~ p / I
            ]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f, p], [I], name=name)
    BondGraph.elements[name] = Element(:I, sys, [sys.p], causality)
    # push!(BondGraph.state_vars, BondGraph.elements[name].p)
end
function add_I!(BondGraph, Φi, params, name; causality=false)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv) p(BondGraph.model.iv)
    D = Differential(BondGraph.model.iv)
    eqns = [
            D(p) ~ e,
            f ~ Φi(p, f, BondGraph.model.iv) # Integral Causality Form
            ]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f, p], [], name=name)
    BondGraph.elements[name] = Element(:I, sys, [sys.p], causality)
    # push!(BondGraph.state_vars, BondGraph.elements[name].p)
end

## Add Memrsistive element
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
function add_Se!(BondGraph, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    @parameters Se(BondGraph.model.iv)
    eqns = [0 ~ e - Se]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], [Se], name=name)
    BondGraph.elements[name] = Element(:Se, sys, [], false)
end

function add_Se!(BondGraph, Se::Number, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [0 ~ e - Se]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], [Se], name=name)
    BondGraph.elements[name] = Element(:Se, sys, [], false)
    push!(BondGraph.inputs, parameters(sys))
end

function add_Se!(BondGraph, Se, params::Vector{}, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [0 ~ e - Se(BondGraph.model.iv, params)]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], params, name=name)
    BondGraph.elements[name] = Element(:Se, sys, [], false)
end

## Add Flow Source
function add_Sf!(BondGraph, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    @parameters Sf(BondGraph.model.iv)
    eqns = [0 ~ f - Sf]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], [Sf], name=name)
    BondGraph.elements[name] = Element(:Se, sys, [], false)
end

function add_Sf!(BondGraph, Sf::Number, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [0 ~ f - Sf]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], [Sf], name=name)
    BondGraph.elements[name] = Element(:Sf, sys, [], false)
    push!(BondGraph.inputs, parameters(sys))
end

function add_Sf!(BondGraph, Sf, params::Vector{}, name)
    @variables e(BondGraph.model.iv) f(BondGraph.model.iv)
    eqns = [0 ~ f - Sf(BondGraph.model.iv, params)]
    sys = ODESystem(eqns, BondGraph.model.iv, [e, f], params, name=name)
    BondGraph.elements[name] = Element(:Sf, sys, [], false)
end

## Add 1-junction to model
function add_1J!(BondGraph, elements::Dict{Symbol,Bool}, name)
    elems = collect(keys(elements))
    eqns = [
            0 ~ sum(x -> BondGraph.elements[x].sys.e * (-1).^(elements[x]), keys(elements)) # Sum of all flows is 0
            ]
    for i ∈ 1:length(elems) - 1
        push!(eqns, BondGraph.elements[elems[i]].sys.f ~ BondGraph.elements[elems[i + 1]].sys.f) # effort equality
    end
    systems = map(x -> BondGraph.elements[x].sys, elems)
    BondGraph.junctions[name] = Junction(:J1, elements, ODESystem(eqns, BondGraph.model.iv), [])
    # return sys
end

## Add 0-junction to model
function add_0J!(BondGraph, elements::Dict{Symbol,Bool}, name)
    eqns = [
            0 ~ sum(x -> BondGraph.elements[x].sys.f * (-1).^(elements[x]), keys(elements)) # Sum of all flows is 0
            ]

    elems = collect(keys(elements))
    for i ∈ 1:length(elems) - 1
        push!(eqns, BondGraph.elements[elems[i]].sys.e ~ BondGraph.elements[elems[i + 1]].sys.e) # effort equality
    end
    # systems = map(x -> BondGraph.elements[x].sys, elems)
    BondGraph.junctions[name] = Junction(:J0, elements, ODESystem(eqns, BondGraph.model.iv, [], [], name=name), [])
end

## Add Transformer to Model
function add_TF!(BondGraph, m, elements::Dict{Symbol,Bool}, name)
    # @parameters m
    elems = collect(keys(elements))
    elem_sys = map(x -> BondGraph.elements[x].sys, elems)
    directions = map(x -> elements[x], elems)
    eqns = [
        0 ~ BondGraph.elements[elems[1]].sys.e * (-1)^directions[1] + m * BondGraph.elements[elems[2]].sys.e * (-1)^directions[2], 
        0 ~ m * BondGraph.elements[elems[1]].sys.f * (-1)^directions[1] + BondGraph.elements[elems[2]].sys.f * (-1)^directions[2]
    ]
    element_sys = map(x -> BondGraph.elements[x].sys, elems)
    sys = ODESystem(eqns, BondGraph.model.iv,  [BondGraph.elements[elems[1]].sys.e, BondGraph.elements[elems[2]].sys.f, BondGraph.elements[elems[1]].sys.f, BondGraph.elements[elems[2]].sys.e], [m], name=name)
    BondGraph.junctions[name] = Junction(:TF, elements, sys, parameters(sys))
end

# Add Gyrator to Model Function add_TF(BondGraph, m, elements::Dict{Symbol, Bool},name)
function add_GY!(BondGraph, r, elements::Dict{Symbol,Bool}, name)
    # @parameters r
    elems = collect(keys(elements))
    directions = map(x -> elements[x], elems)
    eqns = [
        0 ~ BondGraph.elements[elems[1]].sys.e * (-1)^directions[1] + (-1) * r * BondGraph.elements[elems[2]].sys.f * (-1)^directions[2], 
        0 ~ r * BondGraph.elements[elems[1]].sys.f * (-1)^directions[1] + (-1) * BondGraph.elements[elems[2]].sys.e * (-1)^directions[2]
    ]
    element_sys = map(x -> BondGraph.elements[x].sys, elems)
    sys = ODESystem(eqns, BondGraph.model.iv, [BondGraph.elements[elems[1]].sys.e, BondGraph.elements[elems[2]].sys.f, BondGraph.elements[elems[1]].sys.f, BondGraph.elements[elems[2]].sys.e], [r], name=name)
    BondGraph.junctions[name] = Junction(:GY, elements, sys, parameters(sys))
end

# Add MTF - Modulated Transformed Modulus
function add_MTF!(BondGraph, m, states, ps, elements, name)
    # @parameters m
    elems = map(x -> x.first, elements)
    elem_sys = map(x -> BondGraph.elements[x].sys, elems)
    directions = map(x -> x.second, elements)
    eqns = [
        0 ~ BondGraph.elements[elems[1]].sys.e * (-1)^directions[1] + m * BondGraph.elements[elems[2]].sys.e * (-1)^directions[2], 
        0 ~ m * BondGraph.elements[elems[1]].sys.f * (-1)^directions[1] + BondGraph.elements[elems[2]].sys.f * (-1)^directions[2]
    ]
    element_sys = map(x -> BondGraph.elements[x].sys, elems)
    sys = ODESystem(eqns, BondGraph.model.iv,  [BondGraph.elements[elems[1]].sys.e, BondGraph.elements[elems[2]].sys.f, BondGraph.elements[elems[1]].sys.f, BondGraph.elements[elems[2]].sys.e, states...], ps, name=name)
    BondGraph.junctions[name] = Junction(:TF, Dict(elems .=> directions), sys, parameters(sys))
end

## Create ODE System From Bond Graph Construction 
function generate_model!(BondGraph)
    # elements = collect(values(BondGraph.elements))
    # elements = map(x -> x.sys, elements)
    junctions = collect(values(BondGraph.junctions))
    junc_eqns = reduce(vcat, map(x -> equations(x.sys), junctions))
    elem_sys = map(x -> x.sys, collect(values(BondGraph.elements)))
    # state_vars = reduce(vcat,map(x->states(x), elem_sys))
    # elem_sys = reduce(vcat, elem_sys)
    # junc_sys = map(x -> equations(x.sys), collect(values(BondGraph.junctions)))
    # junc_sys = reduce(vcat, junc_sys)
    BondGraph.model = ODESystem(junc_eqns, BondGraph.model.iv, [], [], systems=elem_sys)
end
## Get parameters
function get_parameters!(BondGraph)
    BondGraph.parameters = map(x -> x => 0.0, parameters(BondGraph.model)) |> Dict
    BondGraph.parameters = keys(BondGraph.parameters) .=> values(BondGraph.parameters)
end
## Simplify Bond Graph System 
function simplify_model!(BondGraph)
    BondGraph.model = tearing(structural_simplify(BondGraph.model))
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
<<<<<<< HEAD:BondGraphs/src/BondGraphs.jl
            eqn = expand(simplify(substitute(eqn, sub_dict)))
=======
            eqn = simplify(substitute(eqn, sub_dict))
>>>>>>> fe25e2120125733af19118c21571c6b2ca61d82f:.BondGraphs.jl/src/BondGraphs.jl
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
<<<<<<< HEAD:BondGraphs/src/BondGraphs.jl
        eqn = simplify(substitute(eqn, sub_dict))
=======
        eqn = expand(simplify(substitute(eqn, sub_dict)))
>>>>>>> fe25e2120125733af19118c21571c6b2ca61d82f:.BondGraphs.jl/src/BondGraphs.jl
        push!(res_eqns, eqn)
    end
    return res_eqns
end

function resolve_derivative_causality!(BondGraph)
    # Find the Elements with Derivative Causality aka element.causality=true
    # constiuitive_equations = []
    eqns = equations(BondGraph.model)
    state_vars = reduce(vcat, map(x -> x.state_var, filter(x -> (x.causality == false), collect(values(BondGraph.elements)))))
    D = Differential(BondGraph.model.iv)
    diff_indexes = []
    implicit_eqns = []
    for element ∈ filter(x -> (x.causality == true), collect(values(BondGraph.elements)))
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
    BondGraph.model = ODESystem(eqns, BondGraph.model.iv, states(BondGraph.model), parameters(BondGraph.model))
end

## Start Implementing LTI - Analysis on form ẋ = Ax+Bu, y = Cx +Du
function transfer_function(BondGraph, ps, C, D)
    sts = states(BondGraph.model)
    eqns = equations(BondGraph.model)
    ins = BondGraph.inputs
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

end # module