# using Plots
# using Plots
using OrdinaryDiffEq
using ModelingToolkit: inputs
using SystemBondGraphs
using CairoMakie
using Latexify

# Create a bond graph
@variables t
bg = BondGraph(t)

# Create the Elements in the bond graph
add_R!(bg, :R)
add_C!(bg, :C)
add_I!(bg, :I)
add_Se!(bg, :Se)

# Add the Junction
add_1J!(bg, :J1_1)

# Connect the nodes
add_bond!(bg, :J1_1, :R, :e1)
add_bond!(bg, :J1_1, :C, :e2)
add_bond!(bg, :J1_1, :I, :e3)
add_bond!(bg, :Se, :J1_1, :e4)

# Form the sytem of equations
sys = generate_model(bg)
sys, _ = structural_simplify(sys, (inputs(sys), []))

# Set the inital conditions and parameters
u0 = [sys.C.q => 10.0, sys.I.p => -1.0]
ps = [
    sys.R.R => 1.0,
    sys.C.C => 0.01,
    sys.I.I => 1.0,
    sys.Se.Se => 0.0,
]

# Set the timespan for the simulation
tspan = (0.0, 10.0)

# Create the ODE PRoblem and Solve
prob = ODEProblem(sys, u0, tspan, ps)
sol = solve(prob, Tsit5())
##
latexify(full_equations(sys))
# plot(sol)
##
# set_theme!(theme_latexfonts())
# f = Figure(size = (400,400))
# ax = Axis(f[1,1], xlabel = "Time (t)")
f, ax, p = CairoMakie.plot(sol, axis=(xlabel="Time (t)",))
save("msd.png", f)
# s = sol(0.0:0.01:10.0)
# CairoMakie.lines!(ax, s.t, s[bg[:I].model.p], label = L"p_I")
# CairoMakie.lines!(ax, s.t, s[bg[:C].model.q], label = L"q_C")
# axislegend()
# f
f
