## Add Effort Source
"""
Create a Symbolic/Constant Effort Input. Creates a system with parameters `Se` for the ODESystem.
"""
function add_Se!(BG::AbstractBondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @variables Se(BG.model.iv) [input = true]
    eqns = [0 ~ Se - e]
    model = ODESystem(eqns, BG.model.iv, [e, f, Se], [], name=name)
    type = :Se
    BG.graph[name] = BondGraphNode(name, model, type, [model.Se])
    nothing
end

"""
Create a nonlinear effort input with \$e = S_e(e, f, iv, params)\$ with name
"""
function add_Se!(BG::AbstractBondGraph, Se, params::Vector{}, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ e - Se(e, f, BG.model.iv, params)]
    model = ODESystem(eqns, BG.model.iv, [e, f], params, name=name)
    type = :Se
    BG.graph[name] = BondGraphNode(name, model, type, [])
    nothing
end

## Add Flow Source
"""
Create a Symbolic/Constant Flow Input. Creates a system with parameters `Sf` for the ODESystem.
"""
function add_Sf!(BG::AbstractBondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @variables Sf(BG.model.iv) [input = true]
    eqns = [0 ~ f - Sf]
    model = ODESystem(eqns, BG.model.iv, [e, f, Sf], [], name=name)
    type = :Sf
    BG.graph[name] = BondGraphNode(name, model, type, [model.Sf])
    nothing
end

"""
Create a nonlinear flow input with \$f = S_f(e, f, iv, params)\$ with name
"""
function add_Sf!(BG::AbstractBondGraph, Sf, params::Vector{}, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [0 ~ f - Sf(e, f, BG.model.iv, params)]
    # sys = ODESystem(eqns, BG.model.iv, [e, f], params, name = name)
    model = ODESystem(eqns, BG.model.iv, [e, f], params, name=name)
    type = :Sf
    BG.graph[name] = BondGraphNode(name, model, type, [])
    nothing
end
