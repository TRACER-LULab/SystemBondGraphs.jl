"""

Create a Linear R-Element for the bondgraph. Causality will always be false for an R-element since it does not have a "preferred" causality.

"""
function add_R!(bg, name; causality = false)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv)
    @parameters R
    eqns = [0 ~ R * f - e]
    model = ODESystem(eqns, bg.graph_data.iv, name = name)
    type = :R
    bg[name] = BondGraphNode(model, type)
    nothing
end


"""
Create a Linear C-Element for analysis. Setting Causality to true represents the elements being in derivative causality.
"""
function add_C!(bg, name; causality = false)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv) q(bg.graph_data.iv)
    @parameters C
    D = Differential(bg.graph_data.iv)
    eqns = [D(q) ~ f, 0 ~ q / C - e]
    model = ODESystem(eqns, bg.graph_data.iv, name = name)
    type = :C
    bg[name] = BondGraphNode(model, type)
    nothing
end


"""
Create a Linear I-Element, setting the casuality to true signifies that the element is in derivative causality
"""
function add_I!(bg, name; causality = false)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv) p(bg.graph_data.iv)
    @parameters I
    D = Differential(bg.graph_data.iv)
    eqns = [D(p) ~ e, 0 ~ p / I - f]
    model = ODESystem(eqns, bg.graph_data.iv, name = name)
    type = :I
    bg[name] = BondGraphNode(model, type)

    # add_vertex!(bg.graph)
    # props = Dict(
    #     :type => :I,
    #     :sys => sys,
    #     :causality => causality,
    #     :state_var => [sys.p]
    # )
    # set_prop!(bg.graph, nv(bg.graph), :name, name)
    # set_props!(bg.graph, nv(bg.graph), props)
    nothing
end



"""
Create a Linear M-element with causality being set to false. Derivative causality for M-elements is still under development.
"""
function add_M!(bg, name; causality = false)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv) p(bg.graph_data.iv) q(
        bg.graph_data.iv,
    )
    @parameters M
    D = Differential(bg.graph_data.iv)
    eqns = [D(p) ~ e, D(q) ~ f, 0 ~ M * q - p]
    model = ODESystem(eqns, bg.graph_data.iv, name = name)
    type = :M
    bg[name] = BondGraphNode(model, type)

    # add_vertex!(bg.graph)
    # props = Dict(
    #     :type => :M,
    #     :sys => sys,
    #     :causality => causality,
    #     :state_var => [sys.p, sys.q]
    # )
    # set_prop!(bg.graph, nv(bg.graph), :name, name)
    # set_props!(bg.graph, nv(bg.graph), props)
    nothing
end
