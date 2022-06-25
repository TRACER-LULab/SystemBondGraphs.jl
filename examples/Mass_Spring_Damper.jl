# # Package Imports
using BondGraphs
using DifferentialEquations
using ModelingToolkit
using Plots
# # Create a Bond Graph with independent variable `t`
@parameters t
msd = BondGraph(t, :msd)
# # Add Elements to BondGraph
add_R!(msd, :R_1)
add_C!(msd, :C_1)
add_I!(msd, :I_1)
# # Add Force Input
add_Se!(msd, :Se)
# # Add Junction Connection Elements
add_1J!(msd,
    Dict(
    :R_1 => false,
    :C_1 => false,
    :I_1 => false,
    :Se => true
    ),
    :J1);
# # Traverse Graph and Generate Model
sys = generate_model(msd)
# # Simplify model
model = structural_simplify(sys, simplify = true)
# # Set Initial Condititons and Solve
u0 = [
    msd[:C_1].q => 0.0,
    msd[:I_1].p => 1.0
    ]
ps = [
    msd[:R_1].R => 1.0,
    msd[:C_1].C => 1.0,
    msd[:I_1].I => 1.0,
    msd[:Se].Se => 1.0
    ]
tspan = (0.0, 20.0)
prob = ODAEProblem(model, u0, tspan, ps)
sol = solve(prob, Tsit5())
# # Plot the Results
plot(sol)
