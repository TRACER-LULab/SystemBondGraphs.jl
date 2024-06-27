using OrdinaryDiffEq
using ModelingToolkit: inputs
using SystemBondGraphs
using CairoMakie
using Latexify
##
@variables t
bg = BondGraph(t)
#
add_Se!(bg, :Se)
add_I!(bg, :m)
add_R!(bg, :b1)
add_R!(bg, :b2)
add_C!(bg, :k)
add_Sf!(bg, :vi)

add_1J!(bg, :J1_1)
add_0J!(bg, :J0_1)
add_1J!(bg, :J1_2)
add_0J!(bg, :J0_2)
#
add_bond!(bg, :Se, :J1_1, :e1)
add_bond!(bg, :J1_1, :m, :e2)
add_bond!(bg, :J1_1, :J0_1, :e3)
add_bond!(bg, :J0_1, :b1, :e4)
add_bond!(bg, :J0_1, :J1_2, :e5)
add_bond!(bg, :J1_2, :b2, :e6)
add_bond!(bg, :J1_2, :J0_2, :e7)
add_bond!(bg, :J0_2, :k, :e8)
add_bond!(bg, :vi, :J0_2, :e9)
#
sys = generate_model(bg)
sys, _ = structural_simplify(sys, (inputs(sys), []))
##
ps = [
    sys.m.I => 1.0,
    sys.b1.R => 1.0,
    sys.b2.R => 1.0,
    sys.k.C => 1.0,
    sys.Se.Se => 0.0,
    sys.vi.Sf => 0.0,
]|>Dict
u0 = [
    sys.m.p => 0.0,
    sys.k.q => 10.0,
    sys.b2.e => 5.0,
]|>Dict
tspan = (0.0, 10.0)
prob = ODEProblem(sys, u0, tspan, ps)
sol = solve(prob, Rodas5())
f, ax, p = plot(sol, idxs = [sys.m.p, sys.k.q], axis = (xlabel = "Time (t)",))
axislegend()
f
