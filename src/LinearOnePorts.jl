"""

Create a Linear R-Element for the bondgraph. Causality will always be false for an R-element since it does not have a "preferred" causality.

"""
function add_R!(BG::AbstractBondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters R
    eqns = [0 ~ R * f - e]
    model = ODESystem(eqns, BG.model.iv,name = name)
    type = :R
    BG.graph[name] = BondGraphNode(name, model, type, Num[])
    nothing
end


"""
Create a Linear C-Element for analysis. Setting Causality to true represents the elements being in derivative causality.
"""
function add_C!(BG::AbstractBondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    @parameters C
    D = Differential(BG.model.iv)
    eqns = [
        D(q) ~ f,
        0 ~ q / C - e
    ]
    model = ODESystem(eqns, BG.model.iv,name = name)
    type = :C
    BG.graph[name] = BondGraphNode(name, model, type, [model.q])
    nothing
end


"""
Create a Linear I-Element, setting the casuality to true signifies that the element is in derivative causality
"""
function add_I!(BG::AbstractBondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv)
    @parameters I
    D = Differential(BG.model.iv)
    eqns = [
        D(p) ~ e,
        0 ~ p / I - f
    ]
    model = ODESystem(eqns,BG.model.iv, name = name)
    type = :I
    BG.graph[name] = BondGraphNode(name, model, type, [model.p])

    # add_vertex!(BG.graph)
    # props = Dict(
    #     :type => :I,
    #     :sys => sys,
    #     :causality => causality,
    #     :state_var => [sys.p]
    # )
    # set_prop!(BG.graph, nv(BG.graph), :name, name)
    # set_props!(BG.graph, nv(BG.graph), props)
    nothing
end



"""
Create a Linear M-element with causality being set to false. Derivative causality for M-elements is still under development.
"""
function add_M!(BG::AbstractBondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv) q(BG.model.iv)
    @parameters M
    D = Differential(BG.model.iv)
    eqns = [
        D(p) ~ e,
        D(q) ~ f,
        0 ~ M * q - p
    ]
    model = ODESystem(eqns, BG.model.iv, name = name)
    type = :M
    BG.graph[name] = BondGraphNode(name, model, type, [model.q, model.p])

    # add_vertex!(BG.graph)
    # props = Dict(
    #     :type => :M,
    #     :sys => sys,
    #     :causality => causality,
    #     :state_var => [sys.p, sys.q]
    # )
    # set_prop!(BG.graph, nv(BG.graph), :name, name)
    # set_props!(BG.graph, nv(BG.graph), props)
    nothing
end
