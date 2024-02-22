function add_bond!(bg::AbstractBondGraph, from::Symbol, to::Symbol, name::Symbol)
    @variables e(bg.model.iv) f(bg.model.iv)
    model = ODESystem(Equation[], bg.model.iv, [e,f], [], name=name)
    bg.graph[from, to] = BondGraphEdge(name, model)
    nothing
end
