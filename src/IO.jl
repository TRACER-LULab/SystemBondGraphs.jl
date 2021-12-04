"""
Remove Unecesary fields to only save the graph structure of the bond graph. Returns a graph with only the `:name` and `:type` defined for each node. 
"""
function get_graph(BG::AbstractBondGraph)
    g = BG.graph
    for name in [g.vprops[i][:name] for i ∈ 1:length(keys(g.vprops))]
        rem_prop!(g, g[name, :name],  :causality)
        rem_prop!(g, g[name, :name], :state_var)
        rem_prop!(g, g[name, :name], :sys)
        props(g, g[name, :name]) |> display
    end
    return g
end

"""
Save the BondGraph to a .dot file for later analysis
"""
function savebondgraph(fn::AbstractString, BG::AbstractBondGraph)
    MetaGraphs.savemg(fn, BG.graph)
end

"""
Load and Create a BondGraph from a File and independent variable
"""
function loadbondgraph(fn::AbstractString, iv::Num)
    mg = MetaGraphs.loadmg(fn)
    return BondGraph(mg, iv)
end
