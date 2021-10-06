using BondGraphs
using DifferentialEquations
using ModelingToolkit
using CairoMakie
## Specify Time Parameter
@parameters t
## Create Empty BondGraph
msd = BondGraph(t)
## Add Elements to BondGraph
add_R!(msd, :R_1) 
add_C!(msd, :C_1)
add_I!(msd, :I_1)
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
generate_model!(msd)
## Utilizing Modeling Toolkit to simplify Model
msd.model = structural_simplify(msd.model)
## Normal ModelingToolkit.jl/DifferentialEquations.jl Solving
u0 = [
    msd[:C_1].q => 0.0,
    msd[:I_1].p => -1.0
    ]
ps = [
    msd[:R_1].R => 1.0,
    msd[:C_1].C => 1.0,
    msd[:I_1].I => 1.0,
    msd[:Se].Se => 1.0
    ]
tspan = (0.0, 20.0)
prob = ODAEProblem(msd.model, u0, tspan, ps)
sol = solve(prob, Tsit5())
## Generate Transfer Function
@parameters s
tf = transfer_function(msd, s)
tf[msd[:C_1].q, msd[:Se].Se] |> expand |> simplify