#=
Recreating the work presented in "Analysing and simulating energy-based models in biology using BondGraphTools" by Peter Cudmore, Michael Pan, Peter J. Gawthrop, Edmund J. Crampin 2021
=#
using SystemBondGraphs
using OrdinaryDiffEq
using CairoMakie

# Model Reaction X+A ⇋ Y ⇋ Z ⇋ X+B
@parameters t
bg = ChemBondGraph(t, R=8.314, T=310.0)

## Add Species
add_Ce!(bg, :X)
add_Ce!(bg, :Y)
add_Ce!(bg, :Z)

## Add Zero Junctions
add_0J!(bg, :X0J)
add_0J!(bg, :Y0J)
add_0J!(bg, :Z0J)

## Add Reactions
add_Re!(bg, :XY)
add_Re!(bg, :YZ)
add_Re!(bg, :ZX)

## Add Bonds
add_bond!(bg, :X0J, :X, :e0)
add_bond!(bg, :X0J, :XY, :e1)
add_bond!(bg, :XY, :Y0J, :e2)
add_bond!(bg, :Y0J, :Y, :e3)
add_bond!(bg, :Y0J, :YZ, :e4)
add_bond!(bg, :YZ, :Z0J, :e5)
add_bond!(bg, :Z0J, :Z, :e6)
add_bond!(bg, :Z0J, :ZX, :e7)
add_bond!(bg, :ZX, :X0J, :e8)

## Generate System
sys = generate_model(bg)
sys, _ = structural_simplify(sys, (inputs(sys), []))

## Setup System
## Set Parameters
ps = [
    bg[:X].model.k => 1.0,
    bg[:Y].model.k => 2.0,
    bg[:Z].model.k => 3.0,
    bg[:XY].model.r => 1.0,
    bg[:YZ].model.r => 2.0,
    bg[:ZX].model.r => 3.0,
]

## Set Timespan
tspan = (00.0, 1.0)

## Set Initial Conditions
u0 = [
    bg[:X].model.q => 2.0,
    bg[:Y].model.q => 2.0,
    bg[:Z].model.q => 2.0
]

## Create ODE Problem
prob = ODEProblem(sys, u0, tspan, ps)

## Solve the system
sol = solve(prob, Tsit5())
plot(sol)
