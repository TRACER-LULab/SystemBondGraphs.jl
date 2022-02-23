using BondGraphs
using DifferentialEquations
using ModelingToolkit
## Specify Time Parameter
@parameters t
## Create Empty BondGraph
msd = BondGraph(t)
## Add Elements to BondGraph
add_R!(msd, :R_1) 
add_C!(msd, :C_1)
add_I!(msd, :I_1)
@parameters α ω
# vin(e, f, t, p) = p[1]*sin(t*p[2])
# add_Se!(msd, vin, [α,ω], :Se)
add_Se!(msd, :Se)
## Add Junction Connection Elements 
add_1J!(msd,
    Dict(
    :R_1 => false, 
    :C_1 => false, 
    :I_1 => false, 
    :Se => true
    ),
    :J1);
## Traverse Graph and Generate Model
sys = generate_model(msd)
## Utilizing Modeling Toolkit to simplify Model
model = structural_simplify(sys, simplify = true)
## Normal ModelingToolkit.jl/DifferentialEquations.jl Solving
u0 = [
    msd[:C_1].q => 0.0,
    msd[:I_1].p => 1.0
    ]
ps = [
    msd[:R_1].R => 1.0,
    msd[:C_1].C => 1.0,
    msd[:I_1].I => 1.0,
    msd[:Se].α => 1.0,
    msd[:Se].ω => 1.0,
    ]
tspan = (0.0, 20.0)
prob = ODAEProblem(model, u0, tspan, ps)
sol = solve(prob)