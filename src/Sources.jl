## Add Effort Source
"""
Create a Symbolic/Constant Effort Input. Creates a system with parameters `Se` for the ODESystem.
"""
function add_Se!(BG::AbstractBondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv) 
    @parameters Se(BG.model.iv) [input = true]
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    eqns = [e ~ Se]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Se], name = name)
    props = Dict(
        :type => :Se,
        :sys => sys,
        :causality => false,
        :state_var => [sys.Se]
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

"""
Create a nonlinear effort input with \$e = S_e(e, f, iv, params)\$ with name
"""
function add_Se!(BG::AbstractBondGraph, Se, params::Vector{}, name)
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
        :state_var => []
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

## Add Flow Source
"""
Create a Symbolic/Constant Flow Input. Creates a system with parameters `Sf` for the ODESystem.
"""
function add_Sf!(BG::AbstractBondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters Sf(BG.model.iv) [input = true]
    eqns = [0 ~ f - Sf]
    sys = ODESystem(eqns, BG.model.iv, [e, f], [Sf], name = name)
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    props = Dict(
        :type => :Sf,
        :sys => sys,
        :causality => false,
        :state_var => [sys.Sf]
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

"""
Create a nonlinear flow input with \$f = S_f(e, f, iv, params)\$ with name
"""
function add_Sf!(BG::AbstractBondGraph, Sf, params::Vector{}, name)
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