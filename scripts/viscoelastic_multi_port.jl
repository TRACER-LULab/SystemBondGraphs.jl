using DrWatson
@quickactivate "BondGraphModeling"
# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()
##
using BondGraphs
using ModelingToolkit
using DifferentialEquations
using NLsolve
using NonlinearSolve
## Setup Empty BondGraph
@parameters t
visco = BondGraph(t);
add_Se!(visco, :σ₁)
add_Se!(visco, :σ₂)
add_Se!(visco, :εE₁)
add_Se!(visco, :εE₂)
add_Bond!(visco, :b2)
add_Bond!(visco, :b11)
add_Bond!(visco, :b5)
add_Bond!(visco, :b6)
add_Bond!(visco, :b8)
add_Bond!(visco, :b9)
function ϕi(𝐞, 𝐪, params, BG)
    λ₁, λ₂ = 𝐪
    λ₃ = 1 / λ₁ / λ₂
    I₁ = λ₁^2+λ₂^2+λ₃^2
    μ, Jm = params
    σ1 = (2*λ₁-2*λ₃^2/λ₁)*μ/2/(1-(I₁-3)/Jm)*λ₁
    σ2 = (2*λ₂-2*λ₃^2/λ₂)*μ/2/(1-(I₁-3)/Jm)*λ₂
    return [σ1; σ2]
end 

@parameters μ Jm
elems = [:b5 => true, :b8 => true]
add_C_multiport!(visco, elems, [μ, Jm], :Cα, ϕi = ϕi)

elems = [:b6 => true, :b9 => true]
add_C_multiport!(visco, elems, [μ, Jm], :Cβ, ϕi = ϕi)
add_R!(visco, :R1)
add_R!(visco, :R2)
add_1J!(visco, Dict([
    :σ₁ => true, 
    :εE₁ => true,
    :b5 => false,
    :b2 => false
    ]), :J1_1)
add_1J!(visco, Dict([
    :σ₂ => true, 
    :εE₂ => true,
    :b8 => false,
    :b11 => false
    ]), :J1_2)
add_0J!(visco, Dict([
    :b2 => true,
    :R1 => false,
    :b6 => false
    ]), :J0_1)
add_0J!(visco, Dict([
    :b11 => true,
    :R2 => false,
    :b9 => false
    ]), :J0_2)
##
generate_model!(visco)
visco.model = alias_elimination(visco.model)
##
ps = [
    visco[:Cα].μ  => 18e3,
    visco[:Cβ].μ  => 42e3,
    visco[:Cα].Jm => 110.0,
    visco[:Cβ].Jm => 55.0,
    visco[:R1].R  => 400.0*42e3,
    visco[:R2].R  => 400.0*42e3,
    visco[:σ₁].Se => (ϕi([], [2.0, 2.0], [18e3, 110.0], [])+ϕi([], [2.0, 2.0], [42e3, 55.0], []))[1],
    visco[:σ₂].Se => (ϕi([], [2.0, 2.0], [18e3, 110.0], [])+ϕi([], [2.0, 2.0], [42e3, 55.0], []))[1],
    visco[:εE₁].Se => 0.0,
    visco[:εE₂].Se => 0.0,
]
eqns = equations(visco.model)
eqns = map(x-> substitute(x, Dict(ps))|>simplify, eqns)|>collect

u0 = [
    visco[:Cα].q₁ => 2.0,
    visco[:Cα].q₂ => 2.0,
    visco[:Cβ].q₁ => 1.0,
    visco[:Cβ].q₂ => 1.0,
    visco[:εE₁].e => 0.0,
    visco[:σ₁].e => (ϕi([], [2.0, 2.0], [18e3, 110.0], [])+ϕi([], [2.0, 2.0], [42e3, 55.0], []))[1],
    visco[:b5].e  => (ϕi([], [2.0, 2.0], [18e3, 110.0], [])+ϕi([], [2.0, 2.0], [42e3, 55.0], []))[1],
    visco[:b5].f  => 0.0,
    visco[:εE₂].e => 0.0,
    visco[:σ₂].e => 0.0,
    visco[:b8].e  => 0.0,
    visco[:b8].f  => 0.0,
    visco[:R1].f  => 0.0,
    visco[:b6].f  => 0.0,
    visco[:b6].e  => 0.0,
    visco[:b6].f  => 0.0,
    visco[:R2].f  => 0.0,
    visco[:b9].f  => 0.0,
    visco[:R2].e  => 0.0,
]|>Dict
sys = initialize_system_structure(visco.model)
prob = ODEProblem(visco.model, u0, (0.0, 1.0), ps)
sol = solve(prob)

# ## Generate the model
# nlabels = map(x -> string(get_prop(visco.graph, x, :name)), 1:nv(visco.graph))
# f, ax, p, = graphplot(visco.graph, nlabels = nlabels)
# screen = display(f