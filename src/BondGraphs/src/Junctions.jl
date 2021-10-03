## Add 1-junction to model
function add_1J!(BG::BondGraph, elements::Dict{Symbol,Bool}, name::Symbol)
    elems = collect(keys(elements))
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    for j ∈ keys(elements)
        if elements[j]
            add_edge!(BG.graph, BG.graph[j, :name], BG.graph[name, :name])
        else
            add_edge!(BG.graph, BG.graph[name, :name], BG.graph[j, :name])
        end
    end
    eqns = [
            0 ~ sum(x -> BG[x].e * (-1).^(!elements[x]), keys(elements)) # Sum of all efforts is 0
            ]
    for i ∈ 1:length(elems) - 1
        push!(eqns, BG[elems[i]].f ~ BG[elems[i + 1]].f) # flow equality
    end
    sys = ODESystem(eqns, BG.model.iv, [], [], name = name)
    props = Dict(
                :type => :J1,
                :eqns => eqns,
                :sys  => sys,
                :ps   => []
            )
    set_props!(BG.graph, node_index, props)
    nothing
end

## Add 0-junction to model
function add_0J!(BG::BondGraph, elements::Dict{Symbol,Bool}, name::Symbol)
    elems = collect(keys(elements))
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    for j ∈ keys(elements)
        if elements[j]
            add_edge!(BG.graph, BG.graph[j, :name], BG.graph[name, :name])
        else
            add_edge!(BG.graph, BG.graph[name, :name], BG.graph[j, :name])
        end
    end
    eqns = [
            0 ~ sum(x -> BG[x].f * (-1).^(!elements[x]), keys(elements)) # Sum of all flows is 0
            ]
    for i ∈ 1:length(elems) - 1
        push!(eqns, BG[elems[i]].e ~ BG[elems[i + 1]].e) # effort equality
    end
    sys = ODESystem(eqns, BG.model.iv, [], [], name = name)
    props = Dict(
                :type => :J0,
                :eqns => eqns,
                :sys  => sys,
                :ps   => []
            )
    set_props!(BG.graph, node_index, props)
    nothing
end