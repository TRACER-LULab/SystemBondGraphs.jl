#=
Recreating the work presented in "Analysing and simulating energy-based models in biology using BondGraphTools" by Peter Cudmore, Michael Pan, Peter J. Gawthrop, Edmund J. Crampin 2021
=#
using BondGraphs
using Plots
using ModelingToolkit
using DifferentialEquations
# # Model Reaction X+A ⇋ Y ⇋ Z ⇋ X+B
@parameters t
model = BioBondGraph(t, :chem, R = 8.314, T = 310.0)

# # Build the Model
# ## Add Species
add_Ce!(model, :X)
add_Ce!(model, :Y)
add_Ce!(model, :Z)
# ## Add Building Bonds
for i = 1:6
    add_Bond!(model, Symbol("b" * string(i)))
end
# ## Add Reactions
add_Re!(model, :b1, :b2, :XY)
add_Re!(model, :b3, :b4, :YZ)
add_Re!(model, :b5, :b6, :ZX)
# ## Add 0-Junctions
add_0J!(model,
    Dict(
        :X => false,
        :b6 => true,
        :b1 => false
    ),
    :J0X
)
add_0J!(model,
    Dict(
        :Y => false,
        :b2 => true,
        :b3 => false
    ),
    :J0Y
)
add_0J!(model,
    Dict(
        :Z => false,
        :b4 => true,
        :b5 => false
    ),
    :J0Z
)

# ## Generate Model
sys = generate_model(model)

# ## Setup System
# ### Set Parameters
ps = [
    model[:X].k => 1.0,
    model[:Y].k => 2.0,
    model[:Z].k => 3.0,
    model[:XY].r => 1.0,
    model[:YZ].r => 2.0,
    model[:ZX].r => 3.0,
] |> Dict
# ### Set Timespan
tspan = (00.0, 1.0)
# ### Set Initial Conditions
u0 = [
    model[:X].q => 2.0,
    model[:Y].q => 2.0,
    model[:Z].q => 2.0
] |> Dict
# ### Create ODE Problem
prob = ODEProblem(sys, u0, tspan, ps)

# ## Solve the system
prob = remake(prob, p = [4.0; prob.p[2:end]])
sol = solve(prob, Tsit5())

# ## Plot the solution
plot(sol, ylims = (-2, 4))
