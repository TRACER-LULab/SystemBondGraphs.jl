## Add Transformer to Model
function add_TF!(BG::BondGraph, m, in, out, name)
    # @parameters m
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[in, :name])
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[out, :name])

    eqns = [
        0 ~ BG[in].e - m * BG[out].e, 
        0 ~ m * BG[in].f - BG[out].f
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
function add_GY!(BG, r, in, out, name)
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[in, :name])
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[out, :name])
    eqns = [
        0 ~ BG[in].e - r * BG[out].f, 
        0 ~ r * BG[in].f - BG[out].e
    ]
    sys = ODESystem(eqns, BG.model.iv, name = name)
    add_vertex!(BG.graph)
    props = Dict(
                :type => :GY,
                :eqns => eqns,
                :sys  => sys,
                :ps   => []
            )
    set_props!(BG.graph, node_index, props)  
    nothing
end

## Add MTF - Modulated Transformed Modulus
function add_MTF!(BG::BondGraph, m, in, out, name)
    # @parameters m
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[in, :name])
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[out, :name])

    eqns = [
        0 ~ BG[in].e - m * BG[out].e, 
        0 ~ m * BG[in].f - BG[out].f
    ]
    sys = ODESystem(eqns, BG.model.iv, name = name)
    props = Dict(
                :type => :MTF,
                :eqns => eqns,
                :sys  => sys,
                :ps   => []
            )
    set_props!(BG.graph, node_index, props)
    nothing
end

## Modulated Gyrator
function add_MGY!(BG::BondGraph, r, in, out, name)
    # @parameters m
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[in, :name])
    add_edge!(BG.graph, BG.graph[name, :name], BG.graph[out, :name])

    eqns = [
        0 ~ BG[in].e - r * BG[out].f, 
        0 ~ r * BG[in].f - BG[out].e
    ]
    sys = ODESystem(eqns, BG.model.iv, name = name)
    props = Dict(
                :type => :MTF,
                :eqns => eqns,
                :sys  => sys,
                :ps   => []
            )
    set_props!(BG.graph, node_index, props)
    nothing
end