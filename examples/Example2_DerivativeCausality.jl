using BondGraphs
#
@variables t
bg = BondGraph(t, :bg)
# Add One Ports
add_Se!(bg, :ec)
add_R!(bg, :Rw)
add_I!(bg, :L)
add_C!(bg, :kτ)
add_I!(bg, :J)
add_R!(bg, :bτ)
# Add Two Ports
add_GY!(bg, :T)

# Add Multi Ports
add_1J!(bg, :J1_1)
add_0J!(bg, :J0_1)
add_1J!(bg, :J1_2)

# Connections
add_bond!(bg, :ec, :J1_1, :edge_1)
add_bond!(bg, :J1_1, :L, :edge_2)
add_bond!(bg, :J1_1, :Rw, :edge_3)
add_bond!(bg, :J1_1, :T, :edge_4)
add_bond!(bg, :T, :J0_1, :edge_5)
add_bond!(bg, :J0_1, :kτ, :edge_6)
add_bond!(bg, :J0_1, :J1_2, :edge_7)
add_bond!(bg, :J1_2, :J, :edge_8)
add_bond!(bg, :J1_2, :bτ, :edge_9)
# Generate Model
sys = generate_model(bg)
sys = structural_simplify(sys)
# Simulate System
u0 = [
    bg[:J].p => 10.0,
    bg[:T].e_out => 10.0
]
p = [
    bg[:T].r => 0.5,
    bg[:Rw].R => 0.010,
    bg[:L].I => 1.0,
    bg[:kτ].C => 10.0,
    bg[:J].I => 1.0,
    bg[:bτ].R => 0.1,
    bg[:ec].Se => 20.0,
]
tspan = (0.0, 100.0)
prob = ODAEProblem(sys, u0, tspan, p)
sol = solve(prob, Tsit5())
