using BondGraphs
using ModelingToolkit
using DifferentialEquations

##
@parameters t
mapk = BondGraph(t)

add_Ce!(mapk, :A)
add_Ce!(mapk, :B)
add_Bond!(mapk, :b1)
add_1J!(mapk, Dict(
    :A => true,
    :B => true,
    :b1 => false
    ), :J1)

add_Ce!(mapk, :C)
@parameters re1
add_Re!(mapk, :b1, :C, :Re1)

model = generate_model(mapk)
##
model = structural_simplify(model)

R = 8.3144598
T = 310.0

ps = [
    mapk[:A].k   => 1.0,
    mapk[:B].k   => 0.5,
    mapk[:C].k   => 2.0,
    mapk[:Re1].r => 0.5,
    mapk[:A].R   => R,
    mapk[:B].R   => R,
    mapk[:C].R   => R,
    mapk[:Re1].R => R,
    mapk[:A].T   => T,
    mapk[:B].T   => T,
    mapk[:C].T   => T,
    mapk[:Re1].T => T,
]

u0 = [
    mapk[:A].q => 1.0,
    mapk[:B].q => 2.0,
    mapk[:C].q => 1.0
]
tspan = (0.0, 10.0)
##
prob = ODEProblem(model, u0, tspan, ps)
sol = solve(prob, Tsit5())
# plot(sol)
