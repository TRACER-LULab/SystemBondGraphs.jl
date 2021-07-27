## Add Effort Source
function add_Se!(BG::BondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters Se(BG.model.iv)
    eqns = [0 ~ e - Se]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Se], name = name)
    add_vertex!(BG.graph)
    node_index = length(BG.graph.graph.fadjlist)
    set_prop!(BG.graph, node_index, :name, name)
    set_prop!(BG.graph, node_index, :type, :Se)
    set_prop!(BG.graph, node_index, :sys, sys)
    set_prop!(BG.graph, node_index, :causality, false)
    set_prop!(BG.graph, node_index, :state_var, [])
    nothing
end

function add_Se!(BG::BondGraph, Se::Number, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ e - Se]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Se], name = name)
    push!(BG.inputs, parameters(sys))
    add_vertex!(BG.graph)
    node_index = length(BG.graph.graph.fadjlist)
    set_prop!(BG.graph, node_index, :name, name)
    set_prop!(BG.graph, node_index, :type, :Se)
    set_prop!(BG.graph, node_index, :sys, sys)
    set_prop!(BG.graph, node_index, :causality, false)
    set_prop!(BG.graph, node_index, :state_var, [])
    nothing
end

function add_Se!(BG::BondGraph, Se, params::Vector{}, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ e - Se(BG.model.iv, params)]
    sys = ODESystem(eqns, BG.model.iv, [e, f], params, name = name)
    add_vertex!(BG.graph)
    node_index = length(BG.graph.graph.fadjlist)
    set_prop!(BG.graph, node_index, :name, name)
    set_prop!(BG.graph, node_index, :type, :Se)
    set_prop!(BG.graph, node_index, :sys, sys)
    set_prop!(BG.graph, node_index, :causality, false)
    set_prop!(BG.graph, node_index, :state_var, [])
    nothing
end

## Add Flow Source
function add_Sf!(BG::BondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters Sf(BG.model.iv)
    eqns = [0 ~ f - Sf]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Sf], name = name)
    add_vertex!(BG.graph)
    node_index = length(BG.graph.graph.fadjlist)
    set_prop!(BG.graph, node_index, :name, name)
    set_prop!(BG.graph, node_index, :type, :Sf)
    set_prop!(BG.graph, node_index, :sys, sys)
    set_prop!(BG.graph, node_index, :causality, false)
    set_prop!(BG.graph, node_index, :state_var, [])
    nothing
end

function add_Sf!(BG::BondGraph, Sf::Number, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ f - Sf]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Sf], name = name)
    BG.elements[name] = Element(:Sf, sys, [], false)
    push!(BG.inputs, parameters(sys))
    add_vertex!(BG.graph)
    node_index = length(BG.graph.graph.fadjlist)
    set_prop!(BG.graph, node_index, :name, name)
    set_prop!(BG.graph, node_index, :type, :Sf)
    set_prop!(BG.graph, node_index, :sys, sys)
    set_prop!(BG.graph, node_index, :causality, false)
    set_prop!(BG.graph, node_index, :state_var, [])
    nothing
end

function add_Sf!(BG::BondGraph, Sf, params::Vector{}, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ f - Sf(BG.model.iv, params)]
    sys = ODESystem(eqns, BG.model.iv, [e, f], params, name = name)
    BG.elements[name] = Element(:Sf, sys, [], false)
    add_vertex!(BG.graph)
    node_index = length(BG.graph.graph.fadjlist)
    set_prop!(BG.graph, node_index, :name, name)
    set_prop!(BG.graph, node_index, :type, :Sf)
    set_prop!(BG.graph, node_index, :sys, sys)
    set_prop!(BG.graph, node_index, :causality, false)
    set_prop!(BG.graph, node_index, :state_var, [])
    nothing
end