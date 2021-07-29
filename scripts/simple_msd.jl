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
Makie.inline!(true)
## 
@parameters t
msd = BondGraph(t)
add_R!(msd, :R_1) 
add_C!(msd, :C_1)
add_I!(msd, :I_1)
add_Se!(msd, :Se)
add_1J!(msd,
    Dict(
    :R_1 => false, 
    :C_1 => false, 
    :I_1 => false, 
    :Se => true
    ),
    :J1);
##
generate_model!(msd)
simplify_model!(msd)
## 
u0 = [
    msd[:C_1].q => 0.0,
    msd[:I_1].p => 0.0
    ]
ps = [
    msd[:R_1].R => 1.0,
    msd[:C_1].C => 1.0,
    msd[:I_1].I => 1.0,
    msd[:Se].Se => 1.0
    ]
tspan = (0.0, 10.0)
prob = ODAEProblem(msd.model, u0, tspan, ps)
sol = solve(prob, Tsit5())
