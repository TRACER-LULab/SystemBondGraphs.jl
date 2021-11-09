using BondGraphs
using ModelingToolkit
using DifferentialEquations

##
@parameters t

reaction = BondGraph(t)

add_Ce!(reaction, :A)
add_Ce!(reaction, :B)
add_Re!(reaction, :B, :A, :R13)

model = generate_model(reaction)
model = structural_simplify(model)

@parameters R T r
ps = [
    reaction[:A].k   => 10.0,
    reaction[:A].R   => 1.0,
    reaction[:A].T   => 1.0,
    reaction[:B].k   => 1.0,
    reaction[:B].R   => 1,
    reaction[:B].T   => 1.0,
    r => 1e-1,
    T => 1.0, 
    R => 1.0
]

u0 = [
    reaction[:A].q => 1.0,
    reaction[:B].q => 1.0
]

tspan = (0.0, 5.0)

prob = ODEProblem(model, u0,  tspan, ps)
sol = solve(prob, Rodas4())