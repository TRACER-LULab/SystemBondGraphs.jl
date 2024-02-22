using BondGraphs
using Plots
using DifferentialEquations

##
@parameters t
mapk = BioBondGraph(t, :m, R = 8.214, T = 310.0)

##
add_Ce!(mapk, :A)
add_Ce!(mapk, :B)
add_Ce!(mapk, :C)

##
add_Re!(mapk, :Re1)

##
add_1J!(mapk, :J1)

##
add_bond!(mapk, :A, :J1, :e1)
add_bond!(mapk, :B, :J1, :e2)
add_bond!(mapk, :J1, :Re1, :e3)
add_bond!(mapk, :Re1, :C, :e4)

##
model = generate_model(mapk)
model = structural_simplify(model)

##
ps = [
    mapk[:A].k   => 1.0,
    mapk[:B].k   => 0.5,
    mapk[:C].k   => 2.0,
    mapk[:Re1].r => 0.5,
]

u0 = [
    mapk[:A].q => 1.5,
    mapk[:B].q => 2.0,
    mapk[:C].q => 1.0
]

tspan = (0.0, 10.0)

##
prob = ODEProblem(model, u0, tspan, ps)
sol = solve(prob, Tsit5())
plot(sol)
