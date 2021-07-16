using DrWatson
@quickactivate "BondGraphModeling"
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##
using BondGraphs
using ModelingToolkit
using DifferentialEquations
using LightGraphs
using MetaGraphs
using TikzPictures
using TikzGraphs
##
@variables t
##
msd = BondGraph(t)
##
add_Se!(msd, :Se)
add_R!(msd, :R)
add_C!(msd, :C)
add_I!(msd, :I)
##
add_1J!(msd, Dict(
    :Se => false,
    :R => true,
    :I => true,
    :C => true
    ), :J1_1)
##
name = msd.graph.vprops
name = map(x -> "\$" * string(name[x][:name]) * "\$", 1:length(values(msd.graph.vprops)))
plt = TikzGraphs.plot(msd.graph.graph, name) |> display
TikzGraphs.plot(msd.graph.graph)
##
num_vertices = length(msd.elements) + length(msd.junctions)
g = SimpleGraph(num_vertices)
mg = MetaGraph(g)
set_indexing_prop!(mg, :name)
element_names = collect(keys(msd.elements))
junction_names = collect(keys(msd.junctions))
all_names = [element_names; junction_names]
for i  ∈ eachindex(all_names)
    set_prop!(mg,  i, :name,  all_names[i])
end
for i ∈ keys(msd.junctions)
    for j ∈ keys(msd.junctions[i].elements)
        add_edge!(mg, mg[i, :name], mg[j, :name])
    end
end
plt = TikzGraphs.plot(mg.graph, "\$" .* string.(all_names) .* "\$")
TikzPictures.save(TEX(plotsdir("graph")), plt)