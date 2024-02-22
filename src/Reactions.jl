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
function add_Ce!(BG::BioBondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    @parameters k
    D = Differential(BG.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ BG.R*BG.T*log(k*q)
            ]
    model = ODESystem(
            eqns,
            BG.model.iv,
            [e, f, q],
            [k],
            name = name)
    type = :Ce
    BG.graph[name] = BondGraphNode(name, model, type, [model.q])
    nothing
end

## Biochemical Reaction
# """
# Add a reaction component between species "in" and "out" with parameters R, T, r
# """
# function add_Re!(BG::AbstractBondGraph, in, out, name; causality = false)
#     add_vertex!(BG.graph)
#     set_prop!(BG.graph, nv(BG.graph), :name, name)
#     add_edge!(BG.graph, BG.graph[in, :name], BG.graph[name, :name],)
#     add_edge!(BG.graph, BG.graph[name, :name], BG.graph[out, :name])
#     @parameters r
#     fin = ParentScope(BG[in].f)
#     ein = ParentScope(BG[in].e)
#     fout = ParentScope(BG[out].f)
#     eout = ParentScope(BG[out].e)
#     eqns = [
#         fin ~ fout,
#         fin ~ r*(exp(ein/BG.R/BG.T) - exp(eout/BG.R/BG.T))
#         ]
#     sys  = ODESystem(
#             eqns,
#             BG.model.iv,
#             [],
#             [r],
#             defaults=Dict(r=> 1.0),
#             name = name)
#     props = Dict(
#             :type => :Re,
#             :sys => sys,
#             :causality => causality,
#             :state_var => []
#             )
#     set_props!(BG.graph, nv(BG.graph), props)
#     nothing
# end

"""
Add a reaction component between species "in" and "out" with parameters R, T, r

"""
function add_Re!(bg::BioBondGraph, name)
    @variables e_in(bg.model.iv) f_in(bg.model.iv)
    @variables e_out(bg.model.iv) f_out(bg.model.iv)
    @parameters r

    eqns = [
        f_in ~ f_out,
        f_in ~ r * (exp(e_in / bg.R / bg.T) - exp(e_out / bg.R / bg.T))
    ]

    model = ODESystem(eqns, bg.model.iv, [e_in, e_out, f_in, f_out], [r], name=name)
    type = :Re
    bg.graph[name] = BondGraphNode(name, model, type, Num[])
    nothing
end
