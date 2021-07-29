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
add_C!(mp, :k_1)
add_C!(mp, :k_3; causality = true)
add_Bond!(mp, :b2)
add_Sf!(mp, :F_4)
add_0J!(mp, Dict([
    :F_4 => true,
    :b2 => false,
    :k_3 => false
    ]), :J0)
add_TF!(mp, a / b, :b2, :k_1, :TF)
generate_model!(mp)
equations(mp.model)
resolve_derivative_causality!(mp)

##
simplify_model!(mp)
##
u0 = [
    mp[:k_1].q => 1.0,
    mp[:k_3].q => 0.625,
    mp[:k_3].f => 0.0
    ]
ps = [
    a => 1.0,
    b => 4.0,
    mp[:k_1].C => 1.0,
    mp[:k_3].C => 1 / 4,
    mp[:F_4].Sf => 1.0  
    ]
## 
prob = ODEProblem(mp.model, u0, (0.0, 10.0), ps)
sol = solve(prob)
