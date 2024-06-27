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
add_Se!(bg, :ec)
add_I!(bg, :L)
add_R!(bg, :Rw)
add_GY!(bg, :T)
add_C!(bg, :kτ)
add_R!(bg, :bτ)
add_I!(bg, :J)

# Add the Junction
add_1J!(bg, :J1_1)
add_0J!(bg, :J0_1)
add_1J!(bg, :J1_2)

# Connect the nodes
add_bond!(bg, :ec, :J1_1, :e1)
add_bond!(bg, :J1_1, :Rw, :e2)
add_bond!(bg, :J1_1, :L, :e3)
add_bond!(bg, :J1_1, :T, :e4)
add_bond!(bg, :T, :J0_1, :e5)
add_bond!(bg, :J0_1, :kτ, :e6)
add_bond!(bg, :J0_1, :J1_2, :e7)
add_bond!(bg, :J1_2, :J, :e8)
add_bond!(bg, :J1_2, :bτ, :e9)

# Form the sytem of equations
sys = generate_model(bg)
sys, _ = structural_simplify(sys, (inputs(sys), []))

## Set the inital conditions and parameters
u0 = [
    sys.J.p => 1.0,
    # sys.L.p => 1.0,
    # sys.kτ.q => 1.0,
    ModelingToolkit.D(sys.T.e_out) => 0.0,
]
ps = [
    sys.T.r => 1.0,
    sys.L.I => 1.0,
    sys.Rw.R => 1.0,
    sys.kτ.C => 1.0,
    sys.bτ.R => 1.0,
    sys.J.I => 1.0,
    sys.ec.Se => 0.0,
]

# Set the timespan for the simulation
tspan = (0.0, 10.0)

# Create the ODE PRoblem and Solve
prob = ODEProblem(sys, u0, tspan, ps)
sol = solve(prob, Rodas5())
#
f, ax, p = CairoMakie.plot(sol, idxs = [sys.L.p, sys.kτ.q, sys.J.p], axis=(xlabel="Time (t)",))
save("msd.png", f)
axislegend()
f
