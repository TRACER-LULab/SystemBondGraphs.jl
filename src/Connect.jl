function add_bond!(bg, connection::Pair, name::Symbol)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv)
    model = ODESystem(Equation[], bg.graph_data.iv, [e,f], [], name=name)
    bg[connection.first, connection.second] = BondGraphEdge(name, model)
    nothing
end
