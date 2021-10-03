## Add Effort Source
function add_Se!(BG::BondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters Se(BG.model.iv)
    eqns = [0 ~ e - Se]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Se], name = name)
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    props = Dict(
        :type => :Se,
        :sys => sys,
        :causality => false,
        :state_var => [sys.e]
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

function add_Se!(BG::BondGraph, Se::Number, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ e - Se]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Se], name = name)
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    props = Dict(
        :type => :Se,
        :sys => sys,
        :causality => false,
        :state_var => [sys.e]
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

function add_Se!(BG::BondGraph, Se, params::Vector{}, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ e - Se(e, f, BG.model.iv, params)]
    sys = ODESystem(eqns, BG.model.iv, [e, f], params, name = name)
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    props = Dict(
        :type => :Se,
        :sys => sys,
        :causality => false,
        :state_var => [sys.e]
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

## Add Flow Source
function add_Sf!(BG::BondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters Sf(BG.model.iv)
    eqns = [0 ~ f - Sf]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Sf], name = name)
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    props = Dict(
        :type => :Sf,
        :sys => sys,
        :causality => false,
        :state_var => [sys.f]
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

function add_Sf!(BG::BondGraph, Sf::Number, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ f - Sf]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Sf], name = name)
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    props = Dict(
        :type => :Sf,
        :sys => sys,
        :causality => false,
        :state_var => [sys.f]
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

function add_Sf!(BG::BondGraph, Sf, params::Vector{}, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ f - Sf(e, f, BG.model.iv, params)]
    sys = ODESystem(eqns, BG.model.iv, [e, f], params, name = name)
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    props = Dict(
        :type => :Sf,
        :sys => sys,
        :causality => false,
        :state_var => [sys.f]
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end