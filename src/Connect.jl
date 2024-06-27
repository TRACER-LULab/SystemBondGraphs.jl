function add_bond!(bg, from::Symbol, to::Symbol, name::Symbol)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv)
    model = ODESystem(Equation[], bg.graph_data.iv, [e,f], [], name=name)
    bg[from, to] = BondGraphEdge(name, model)
    nothing
end
