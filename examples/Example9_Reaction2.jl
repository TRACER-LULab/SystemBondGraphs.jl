using BondGraphs
using OrdinaryDiffEq
using CairoMakie
# Plot settings
set_theme!(theme_latexfonts())
f = Figure(size=(350, 350))
ax = Axis(f[1, 1], xlabel="Time (t)", ylabel="Concentration (q)")
# Create the chemical bond graph
@parameters t
mapk = ChemBondGraph(t, R = 1.0, T = 1.0)

# Add Species
add_Ce!(mapk, :A)
add_Ce!(mapk, :B)
add_Ce!(mapk, :C)

# Add Reaction
add_Re!(mapk, :Re)

# Combine the species and the reaction
add_1J!(mapk, :J1)

# Add Bonds
add_bond!(mapk, :A, :J1, :e1)
add_bond!(mapk, :B, :J1, :e2)
add_bond!(mapk, :J1, :Re, :e3)
add_bond!(mapk, :Re, :C, :e4)

# Generate model and simplify
model = generate_model(mapk)
model = structural_simplify(model)

# Set parameters
ps = [
    mapk[:A].model.k   => 1.0,
    mapk[:B].model.k   => 1.0,
    mapk[:C].model.k   => 1.0,
    mapk[:Re].model.r => 1.0,
]
# Set inital conditions
u0 = [
    mapk[:A].model.q => 0.0,
    mapk[:B].model.q => 1.0,
    mapk[:C].model.q => 1.0
]
# Set time span
tspan = (0.0, 2.0)

# Solve the ODE
prob = ODEProblem(model, u0, tspan, ps)
sol = solve(prob, Tsit5())
# Plot the results
s = sol(0.0:0.01:2.0)
CairoMakie.lines!(ax, s.t, s[mapk[:A].model.q], label=L"q_A")
CairoMakie.lines!(ax, s.t, s[mapk[:B].model.q], label=L"q_B")
CairoMakie.lines!(ax, s.t, s[mapk[:C].model.q], label=L"q_C")
axislegend(position=:rb)
f
