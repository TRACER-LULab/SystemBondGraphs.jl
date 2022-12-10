using BondGraphs
using ModelingToolkit
using DifferentialEquations
# # Create the bond graph
@parameters t
reaction = BioBondGraph(t, :re1)
# # Add Species and reaction
add_Ce!(reaction, :A)
add_Ce!(reaction, :B)
add_Re!(reaction, :B, :A, :R13)
# # Generate Mode
model = generate_model(reaction)
model = structural_simplify(model)
# Set Parameters
@parameters R T r
ps = [
    reaction[:A].k   => 10.0,
    reaction[:B].k   => 1.0,
    reaction[:R13].r => 0.1,
]
u0 = [
    reaction[:A].q => 10.0,
    reaction[:B].q => 1.0
]
tspan = (0.0, 5.0)
# # Solve and Plot
prob = ODEProblem(model, u0,  tspan, ps)
sol = solve(prob)
plot(sol)
