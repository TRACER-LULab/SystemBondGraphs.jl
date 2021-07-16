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
## Plot Underlying Graph
name = msd.graph.vprops
name = map(x -> "\$" * string(name[x][:name]) * "\$", 1:length(values(msd.graph.vprops)))
plt = TikzGraphs.plot(msd.graph.graph, name) |> display
