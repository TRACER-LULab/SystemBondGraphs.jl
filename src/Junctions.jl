function add_1J!(bg, name::Symbol)
    model = ODESystem(Equation[],bg.graph_data.iv, name = name)
    type = :OneJunction
    bg[name] = BondGraphNode(model, type)
    nothing
end

function add_0J!(bg, name::Symbol)
    model = ODESystem(Equation[], bg.graph_data.iv, name=name)
    type = :ZeroJunction
    bg[name] = BondGraphNode(model, type)
    nothing
end
