using SystemBondGraphs
using OrdinaryDiffEq
using CairoMakie
# Create the bond graph
@parameters t, R, T
bg = ChemBondGraph(t, R=1.0, T=300.0)
# Add Species and reaction
add_Ce!(bg, :A)
add_Ce!(bg, :B)
add_Re!(bg, :R13)
# Add Bonds
add_bond!(bg, :B, :R13, :e1)
add_bond!(bg, :R13, :A, :e2)

# Generate Mode
sys = generate_model(bg)
sys, _ = structural_simplify(sys, (inputs(sys), []))
# Set Parameters
ps = [
    bg[:A].model.k => 0.10,
    bg[:B].model.k => 1.0,
    bg[:R13].model.r => 1.0,
]
u0 = [bg[:A].model.q => 1.0, bg[:B].model.q => 1.0]
tspan = (0.0, 10.0)

# # Solve and Plot
prob = ODEProblem(sys, u0, tspan, ps)
sol = solve(prob, Tsit5())
f, ax, p = plot(sol)
axislegend(position = :rc)
f
