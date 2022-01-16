## Add 1-junction to model
"""
Add a 1-Junction to the bond graph, `BG`. The elements are specified with a dict of names and directions. `true` represents the element going "in" to the 1-junction. `false` represents the element going "out" from the 1-junction
"""
function add_1J!(BG::AbstractBondGraph, elements::Dict{Symbol,Bool}, name::Symbol)
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
    sys = ODESystem(Equation[], BG.model.iv, [], [], name = name)
    props = Dict(
                :type => :J1,
                :sys  => sys,
            )
    set_props!(BG.graph, node_index, props)
    nothing
end

## Add 0-junction to model
"""
Add a 0-Junction to the bond graph, `BG`. The elements are specified with a dict of names and directions. `true` represents the element going "in" to the 0-junction. `false` represents the element going "out" from the 0-junction
"""
function add_0J!(BG::AbstractBondGraph, elements::Dict{Symbol,Bool}, name::Symbol)
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
    sys = ODESystem(Equation[], BG.model.iv, [], [], name = name)
    props = Dict(
                :type => :J0,
                :sys  => sys,
            )
    set_props!(BG.graph, node_index, props)
    nothing
end