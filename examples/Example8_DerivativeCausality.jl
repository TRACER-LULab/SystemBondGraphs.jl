using SystemBondGraphs
using OrdinaryDiffEq
using CairoMakie
#
@variables t
bg = BondGraph(t)
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
sys, _ = structural_simplify(sys, (inputs(sys), []))
# Simulate System
e_out_t = ModelingToolkit.D(bg[:T].model.e_out)
u0 = [bg[:J].model.p => 10.0, bg[:kτ].model.e => 10.0, e_out_t => 0.0] |> Dict
p =
    [
        bg[:T].model.r => 0.5,
        bg[:Rw].model.R => 0.010,
        bg[:L].model.I => 1.0,
        bg[:kτ].model.C => 10.0,
        bg[:J].model.I => 1.0,
        bg[:bτ].model.R => 0.1,
        bg[:ec].model.Se => 20.0,
    ] |> Dict
tspan = (0.0, 100.0)
prob = ODEProblem(sys, u0, tspan, p)
sol = solve(prob, Rodas5())
CairoMakie.plot(sol)
