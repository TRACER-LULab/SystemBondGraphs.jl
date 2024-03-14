## Add Effort Source
"""
Create a Symbolic/Constant Effort Input. Creates a system with parameters `Se` for the ODESystem.
"""
function add_Se!(bg, name)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv)
    @variables Se(bg.graph_data.iv) [input = true]
    eqns = [0 ~ Se - e]
    model = ODESystem(eqns, bg.graph_data.iv, [e, f, Se], [], name=name)
    type = :Se
    bg[name] = BondGraphNode(model, type)
    nothing
end

"""
Create a nonlinear effort input with \$e = S_e(e, f, iv, params)\$ with name
"""
function add_Se!(bg, Se, params::Vector{}, name)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv)
    eqns = [0 ~ e - Se(e, f, bg.graph_data.iv, params)]
    model = ODESystem(eqns, bg.graph_data.iv, [e, f], params, name=name)
    type = :Se
    bg[name] = BondGraphNode(model, type)
    nothing
end

## Add Flow Source
"""
Create a Symbolic/Constant Flow Input. Creates a system with parameters `Sf` for the ODESystem.
"""
function add_Sf!(bg, name)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv)
    @variables Sf(bg.graph_data.iv) [input = true]
    eqns = [0 ~ f - Sf]
    model = ODESystem(eqns, bg.graph_data.iv, [e, f, Sf], [], name=name)
    type = :Sf
    bg[name] = BondGraphNode(model, type)
    nothing
end

"""
Create a nonlinear flow input with \$f = S_f(e, f, iv, params)\$ with name
"""
function add_Sf!(bg, Sf, params::Vector{}, name)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv)
    eqns = [0 ~ f - Sf(e, f, bg.graph_data.iv, params)]
    # sys = ODESystem(eqns, bg.graph_data.iv, [e, f], params, name = name)
    model = ODESystem(eqns, bg.graph_data.iv, [e, f], params, name=name)
    type = :Sf
    bg[name] = BondGraphNode(model, type)
    nothing
end
