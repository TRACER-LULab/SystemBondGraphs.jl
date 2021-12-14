## Add Chemostats - BioModeling
export add_Chemostat!,
        add_Flowstat!,
        add_Ce!,
        add_Re!
"""
Create a linear chemostat
"""
add_Chemostat!(BG::AbstractBondGraph, name) = add_Se!(BG, name)
"""
Create a nonliner chemostat
"""
add_Chemostat!(BG::AbstractBondGraph, Se, params::Vector{}, name) = add_Se!(BG::AbstractBondGraph, Se, params::Vector{}, name)

## Add a Flowstat
"""
Create a linear Flowstat
"""
add_Flowstat!(BG::AbstractBondGraph, name) = add_Sf!(BG, name)
"""
Create a nonlinear Flowstat
"""
add_Flowstat!(BG::AbstractBondGraph, Se, params::Vector{}, name) = add_Sf!(BG::AbstractBondGraph, Se, params::Vector{}, name)

## Concentration Element?
"""
Create a "Linear" Concentration-Element for analysis. Setting Causality to true represents the elements being in derivative causality. 
"""
function add_Ce!(BG::AbstractBondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    @parameters k
    D = Differential(BG.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ BG.R*BG.T*log(k*q)
            ]
    sys = ODESystem(
            eqns, 
            BG.model.iv, 
            [e, f, q], 
            [k], 
            defaults = Dict(q=>0.0),
            name = name)
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
function add_Re!(BG::AbstractBondGraph, in, out, name; causality = false)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    add_edge!(BG.graph, BG.graph[in, :name], BG.graph[name, :name],)
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[out, :name])    
    @parameters r
    fin = ParentScope(BG[in].f)
    ein = ParentScope(BG[in].e)
    fout = ParentScope(BG[out].f)
    eout = ParentScope(BG[out].e)
    eqns = [
        fin ~ fout,
        fin ~ r*(exp(ein/BG.R/BG.T) - exp(eout/BG.R/BG.T))
        ]
    sys  = ODESystem(
            eqns, 
            BG.model.iv, 
            [], 
            [r], 
            defaults=Dict(r=> 1.0),
            name = name)
    props = Dict(
            :type => :Re,
            :sys => sys,
            :causality => causality,
            :state_var => []
            )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end