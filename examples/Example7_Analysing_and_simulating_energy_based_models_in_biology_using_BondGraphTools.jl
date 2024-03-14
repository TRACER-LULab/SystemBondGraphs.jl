#=
Recreating the work presented in "Analysing and simulating energy-based models in biology using BondGraphTools" by Peter Cudmore, Michael Pan, Peter J. Gawthrop, Edmund J. Crampin 2021
=#
using BondGraphs
using OrdinaryDiffEq
using Plots

# Model Reaction X+A ⇋ Y ⇋ Z ⇋ X+B
@parameters t
model = BioBondGraph(t, R = 8.314, T = 310.0)

## Add Species
add_Ce!(model, :X)
add_Ce!(model, :Y)
add_Ce!(model, :Z)

## Add Zero Junctions
add_0J!(model, :X0J)
add_0J!(model, :Y0J)
add_0J!(model, :Z0J)

## Add Reactions
add_Re!(model, :XY)
add_Re!(model, :YZ)
add_Re!(model, :ZX)

## Add Bonds
add_bond!(model, :X0J, :X, :e0)
add_bond!(model, :X0J, :XY, :e1)
add_bond!(model, :XY, :Y0J, :e2)
add_bond!(model, :Y0J, :Y, :e3)
add_bond!(model, :Y0J, :YZ, :e4)
add_bond!(model, :YZ, :Z0J, :e5)
add_bond!(model, :Z0J, :Z, :e6)
add_bond!(model, :Z0J, :ZX, :e7)
add_bond!(model, :ZX, :X0J, :e8)

## Generate System
sys = generate_model(model)
sys = structural_simplify(sys)

## Setup System
## Set Parameters
ps = [
    model[:X].model.k => 1.0,
    model[:Y].model.k => 2.0,
    model[:Z].model.k => 3.0,
    model[:XY].model.r => 1.0,
    model[:YZ].model.r => 2.0,
    model[:ZX].model.r => 3.0,
]

## Set Timespan
tspan = (00.0, 1.0)

## Set Initial Conditions
u0 = [
    model[:X].model.q => 2.0,
    model[:Y].model.q => 2.0,
    model[:Z].model.q => 2.0
]

## Create ODE Problem
prob = ODEProblem(sys, u0, tspan, ps)

## Solve the system
sol = solve(prob, Tsit5())
plot(sol)
