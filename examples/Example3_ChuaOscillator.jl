# # Nonlinear Chaotic Chua Oscillator
using BondGraphs
using DifferentialEquations
# Make the BondGraph
@variables t
bg = BondGraph(t, :bg)
# Add Nonlinear Resistor Equation
@parameters Ga Gb E
Φr(e, f, t, p) = p[2] * e + 1 / 2 * (p[1] - p[2]) * (abs(e + p[3]) - abs(e - p[3]))
add_R!(bg, :f => Φr, [Ga, Gb, E], :Gn)
# Add Linear Elements
add_R!(bg, :G)
add_R!(bg, :R)
add_C!(bg, :C1)
add_C!(bg, :C2)
add_I!(bg, :L)
# Add Junctions
add_1J!(bg, :J11)
add_1J!(bg, :J12)
add_0J!(bg, :J01)
add_1J!(bg, :J13)
add_0J!(bg, :J02)
# Add Connections
add_bond!(bg, :J11, :C2, :edge_1)
add_bond!(bg, :J11, :J01, :edge_2)
add_bond!(bg, :J12, :G, :edge_3)
add_bond!(bg, :J12, :J01, :edge_4)
add_bond!(bg, :J01, :J13, :edge_5)
add_bond!(bg, :J13, :L, :edge_6)
add_bond!(bg, :J13, :R, :edge_7)
add_bond!(bg, :J13, :J02, :edge_8)
add_bond!(bg, :J02, :C1, :edge_9)
add_bond!(bg, :J02, :Gn, :edge_10)
# Generate System
model = generate_model(bg)
model = structural_simplify(model)
u0 = [
    bg[:C1].q => 0.0,
    bg[:C2].q => -0.1,
    bg[:L].p => 0.0001
]
ps = [
    bg[:R].R => 12.5e-3,
    bg[:G].R => 0.565,
    bg[:C1].C => 10.0,
    bg[:C2].C => 100.0,
    bg[:L].I => 18.0,
    bg[:Gn].E => -2.0,
    bg[:Gn].Ga => 0.757576,
    bg[:Gn].Gb => 0.409091
]
tspan = (0.0, 1e3)
prob = ODEProblem(model, u0, tspan, ps)
sol = solve(prob, Tsit5())
