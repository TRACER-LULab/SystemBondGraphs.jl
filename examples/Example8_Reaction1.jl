using BondGraphs
using OrdinaryDiffEq
using Plots
# Create the bond graph
@parameters t
reaction = ChemBondGraph(t)
# Add Species and reaction
add_Ce!(reaction, :A)
add_Ce!(reaction, :B)
add_Re!(reaction, :R13)
# Add Bonds
add_bond!(reaction, :B, :R13, :e1)
add_bond!(reaction, :R13, :A, :e2)

# Generate Mode
model = generate_model(reaction)
model = structural_simplify(model)
# Set Parameters
ps = [
    reaction[:A].model.k   => 0.10,
    reaction[:B].model.k   => 1.0,
    reaction[:R13].model.r => 1.0,
]
u0 = [
    reaction[:A].model.q => 1.0,
    reaction[:B].model.q => 1.0
]
tspan = (0.0, 10.0)

# # Solve and Plot
prob = ODEProblem(model, u0,  tspan, ps)
sol = solve(prob, Tsit5())
Plots.plot(sol)
