## Swap new node for old node and delete old node. 
function swap(BG, old, new)
    old = S(old)
    new = S(new)
    in_nodes = inneighbors(BG.graph, BG.graph[old, :name])
    out_nodes = outneighbors(BG.graph, BG.graph[old, :name])
    if !isempty(in_nodes)
        for i ∈ in_nodes
            add_edge!(BG.graph, i, BG.graph[new, :name])
        end
    end
    if !isempty(out_nodes)
        for i ∈ out_nodes
            add_edge!(BG.graph, BG.graph[new, :name], i)
        end
    end
    rem_vertex!(BG.graph, BG.graph[old, :name])
end
## Apply defaults from "old" to "new" element of the same type
function swap_defaults(BG, name, old, new)
    for key in collect(keys(BG[S(name * old)].defaults))
        BG[S(name * new)].defaults[key] = BG[S(name * old)].defaults[key]
    end
end
## Swap element for a bond element
function swap_bond(BG, elem)
    old_name = S(elem)
    new_name = S(elem * "_bond")
    add_Bond!(BG, new_name)
    swap(BG, elem, string(new_name))
    set_prop!(BG.graph, BG.graph[new_name, :name], :name, old_name)
end
