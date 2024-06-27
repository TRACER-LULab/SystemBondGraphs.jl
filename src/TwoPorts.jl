## Add Transformer to Model\
"""
Add a linear transformer with modulus `m`, in element name `in`, out element name `out`, and named `name`.
"""
function add_TF!(bg, name)
    @variables e_in(bg.graph_data.iv) f_in(bg.graph_data.iv)
    @variables e_out(bg.graph_data.iv) f_out(bg.graph_data.iv)
    @parameters m
    eqns = [
        0 ~ e_in - m * e_out,
        0 ~ m * f_in - f_out
    ]
    model = ODESystem(eqns, bg.graph_data.iv, [e_in, e_out, f_in, f_out], [m], name=name)
    type = :TF
    bg[name] = BondGraphNode(model, type)
    nothing
end

## Add Gyrator to Model Function add_TF(bg, m, elements::Dict{Symbol, Bool},name)
"""
Add a linear gyrator with modulus `r`, in element name `in`, out element name `out`, and named `name`.
"""
function add_GY!(bg, name)

    @variables e_in(bg.graph_data.iv) f_in(bg.graph_data.iv)
    @variables e_out(bg.graph_data.iv) f_out(bg.graph_data.iv)
    @parameters r [description="Modulus of the gyrator"]

    eqns = [
        0 ~ e_in - r * f_out,
        0 ~ r * f_in - e_out
    ]
    model = ODESystem(eqns, bg.graph_data.iv, [e_in, e_out, f_in, f_out], [r], name=name)
    type = :GY
    bg[name] = BondGraphNode(model, type)
    nothing
end

## Add MTF - Modulated Transformed Modulus
"""
Add a modulated transformer with modulus `m`, in element name `in`, out element name `out`, and named `name`.
"""
function add_MTF!(bg, m, name)
    # @parameters m
    # add_vertex!(bg.graph)
    # node_index = nv(bg.graph)
    # set_prop!(bg.graph, node_index, :name, name)
    # add_edge!(bg.graph, bg.graph[name, :name], bg.graph[in, :name])
    # add_edge!(bg.graph, bg.graph[name, :name], bg.graph[out, :name])
    @variables e_in(bg.graph_data.iv) f_in(bg.graph_data.iv)
    @variables e_out(bg.graph_data.iv) f_out(bg.graph_data.iv)

    eqns = [
        0 ~ e_in - ParentScope(m) * e_out,
        0 ~ ParentScope(m) * f_in - f_out
    ]
    model = ODESystem(eqns, bg.graph_data.iv, name=name)
    type = :MTF
    bg[name] = BondGraphNode(model, type)

    # props = Dict(
    #             :type => :MTF,
    #             :eqns => eqns,
    #             :sys  => sys,
    #             :ps   => []
    #         )
    # set_props!(bg.graph, node_index, props)
    nothing
end

## Modulated Gyrator
"""
Add a modulated gyrator with modulus `r`, in element name `in`, out element name `out`, and named `name`.
"""
function add_MGY!(bg, r, name)
    @variables e_in(bg.graph_data.iv) f_in(bg.graph_data.iv)
    @variables e_out(bg.graph_data.iv) f_out(bg.graph_data.iv)

    eqns = [
        0 ~ e_in - r * f_out,
        0 ~ r * f_in - e_out
    ]

    model = ODESystem(eqns, bg.graph_data.iv, name=name)
    type = :MGY
    bg[name] = BondGraphNode(model, type)

    nothing
end
