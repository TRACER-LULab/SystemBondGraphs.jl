function add_1J!(bg::AbstractBondGraph, name::Symbol)
    model = ODESystem(Equation[],bg.model.iv, name = name)
    type = :OneJunction
    bg.graph[name] = BondGraphNode(name, model, type, Num[])
    nothing
end

function add_0J!(bg::AbstractBondGraph, name::Symbol)
    model = ODESystem(Equation[], bg.model.iv, name=name)
    type = :ZeroJunction
    bg.graph[name] = BondGraphNode(name, model, type, Num[])
    nothing
end
