## Add Transformer to Model\
"""
Add a linear transformer with modulus `m`, in element name `in`, out element name `out`, and named `name`.
"""
function add_TF!(BG::AbstractBondGraph, name)
    # @parameters m
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[in, :name])
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[out, :name])
    @parameters m
    eqns = [
        0 ~ ParentScope(BG[in].e) - m * ParentScope(BG[out].e),
        0 ~ m * ParentScope(BG[in].f) - ParentScope(BG[out].f)
    ]

    sys = ODESystem(eqns, BG.model.iv, name = name)
    props = Dict(
                :type => :TF,
                :eqns => eqns,
                :sys  => sys,
                :ps   => []
            )
    set_props!(BG.graph, node_index, props)
    nothing
end

## Add Gyrator to Model Function add_TF(BG, m, elements::Dict{Symbol, Bool},name)
"""
Add a linear gyrator with modulus `r`, in element name `in`, out element name `out`, and named `name`.
"""
function add_GY!(bg::AbstractBondGraph, name)

    @variables e_in(bg.model.iv) f_in(bg.model.iv)
    @variables e_out(bg.model.iv) f_out(bg.model.iv)
    @parameters r

    eqns = [
        0 ~ e_in - r * f_out,
        0 ~ r * f_in - e_out
    ]
    model = ODESystem(eqns, bg.model.iv, [e_in, e_out, f_in, f_out], [r], name=name)
    type = :GY
    bg.graph[name] = BondGraphNode(name, model, type, Num[])
    nothing
end

## Add MTF - Modulated Transformed Modulus
"""
Add a modulated transformer with modulus `m`, in element name `in`, out element name `out`, and named `name`.
"""
function add_MTF!(bg::AbstractBondGraph, m, name)
    # @parameters m
    # add_vertex!(BG.graph)
    # node_index = nv(BG.graph)
    # set_prop!(BG.graph, node_index, :name, name)
    # add_edge!(BG.graph, BG.graph[name, :name], BG.graph[in, :name])
    # add_edge!(BG.graph, BG.graph[name, :name], BG.graph[out, :name])
    @variables e_in(bg.model.iv) f_in(bg.model.iv)
    @variables e_out(bg.model.iv) f_out(bg.model.iv)

    eqns = [
        0 ~ e_in - ParentScope(m) * e_out,
        0 ~ ParentScope(m) * f_in - f_out
    ]
    model = ODESystem(eqns, bg.model.iv, name = name)
    type = :MTF
    bg.graph[name] = BondGraphNode(name, model, type, Num[])

    # props = Dict(
    #             :type => :MTF,
    #             :eqns => eqns,
    #             :sys  => sys,
    #             :ps   => []
    #         )
    # set_props!(BG.graph, node_index, props)
    nothing
end

## Modulated Gyrator
"""
Add a modulated gyrator with modulus `r`, in element name `in`, out element name `out`, and named `name`.
"""
function add_MGY!(BG::AbstractBondGraph, r, in, out, name)
    # @parameters m
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[in, :name])
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[out, :name])

    eqns = [
        0 ~ ParentScope(BG[in].e) - r * ParentScope(BG[out].f),
        0 ~ r * ParentScope(BG[in].f) - ParentScope(BG[out].e)
    ]
    sys = ODESystem(eqns, BG.model.iv, name = name)
    props = Dict(
        :type => :MTF,
        :eqns => eqns,
        :sys => sys,
        :ps => []
    )
    set_props!(BG.graph, node_index, props)
    nothing
end
