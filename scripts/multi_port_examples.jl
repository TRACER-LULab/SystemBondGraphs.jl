using DrWatson
@quickactivate "BondGraphModeling"
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##  
using BondGraphs
using DifferentialEquations
using ModelingToolkit
## 
@parameters t a b
mp = BondGraph(t)
## Fig 7.8a
add_C!(mp, :k_1)
add_C!(mp, :k_2; causality = true)
add_Bond!(mp, :b2)
add_Se!(mp, :F_4)
##
add_0J!(mp, Dict([
    :F_4 => true,
    :b2 => false,
    :k_2 => false
    ]), :J0)

add_TF!(mp, a / b, :b2, :k_1, :TF)
##
generate_model!(mp)
resolve_derivative_causality!(mp)
simplify_model!(mp)
##
