using DrWatson
@quickactivate "BondGraphModeling"
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##  
using BondGraphs
using DifferentialEquations
using ModelingToolkit
using LightGraphs
using MetaGraphs
using GLMakie
using GraphMakie
using LaTeXStrings
## 
@parameters t
msd = BondGraph(t)
##
add_R!(msd, :R_1) 
add_C!(msd, :C_1)
add_I!(msd, :I_1)
add_Se!(msd, :Se)
add_1J!(msd, Dict([
    :R_1 => true, 
    :C_1 => true, 
    :I_1 => true, 
    :Se => false]),
    :J1)
##
generate_model!(msd)
##
names = [latexstring(msd.graph.vprops[i][:name]) for i in 1:length(msd.graph.vprops)]
f, ax, p = graphplot(g, nlabels=names)
offsets = 0.02 * (p[:node_positions][] .- p[:node_positions][][1])
p.nlabels_offset[] = offsets
autolimits!(ax)