## Add Chemostats - BioModeling
export add_Chemostat!,
        add_Flowstat!,
        add_Ce!,
        add_Re!
"""
Create a linear chemostat
"""
add_Chemostat!(BG::BondGraph, name) = add_Se!(BG, name)
"""
Create a nonliner chemostat
"""
add_Chemostat!(BG::BondGraph, Se, params::Vector{}, name) = add_Se!(BG::BondGraph, Se, params::Vector{}, name)

## Add a Flowstat
"""
Create a linear Flowstat
"""
add_Flowstat!(BG::BondGraph, name) = add_Sf!(BG, name)
"""
Create a nonlinear Flowstat
"""
add_Flowstat!(BG::BondGraph, Se, params::Vector{}, name) = add_Sf!(BG::BondGraph, Se, params::Vector{}, name)

## Concentration Element?
"""
Create a "Linear" Concentration-Element for analysis. Setting Causality to true represents the elements being in derivative causality. 
"""
function add_Ce!(BG::BondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    @parameters R T k
    D = Differential(BG.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ R*T*log(k*q)
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, q], [R, T, k], name = name)
    add_vertex!(BG.graph)
    props = Dict(
            :type => :Ce,
            :sys => sys,
            :causality => causality,
            :state_var => [sys.q]
            )
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

## Biochemical Reaction
"""
Add a reaction component between species "in" and "out" with parameters R, T, r
"""
function add_Re!(BG::BondGraph, in, out, name; causality = false)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    @parameters R T r
    eqns = [
        0 ~ BG[in].f + BG[out].f,
        0 ~ BG[in].f + r*(exp(BG[in].e/R/T) - exp(BG[out].e/R/T))
        ]
    sys  = ODESystem(eqns, BG.model.iv, [], [R, T, r], name = name)
    props = Dict(
            :type => :Re,
            :sys => sys,
            :causality => causality,
            :state_var => []
            )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end