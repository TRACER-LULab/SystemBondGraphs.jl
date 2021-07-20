using DrWatson
@quickactivate "BondGraphModeling"
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##
using BondGraphs
using ModelingToolkit
using DifferentialEquations
using LightGraphs
using MetaGraphs
using TikzPictures
using TikzGraphs
##
@variables t
##
msd = BondGraph(t)
##
add_Se!(msd, :Se)
add_R!(msd, :R)
add_C!(msd, :C)
add_I!(msd, :I)
add_M!(msd, :M)
##
add_1J!(msd, Dict(
    :Se => false,
    :R => true,
    :I => true,
    :C => true,
    :M => true
    ), :J1_1)
## 
generate_model!(msd)
structural_simplify(msd.model)

tg = BondGraph(t)

function s(𝐞, 𝐪, params)
    λ₁, λ₂ = 𝐪
    μ = params
    s1 = μ / 2 * (2 * λ₁ - 2 * λ₁^(-3) * λ₂^(-2))
    s2 = μ / 2 * (2 * λ₂ - 2 * λ₁^(-2) * λ₂^(-3))
    return [s1; s2]
end 

add_Bond!(tg, :b1)
add_Bond!(tg, :b2)

elems = [:b1 => false, :b2 => false]

add_C_multiport!(tg, elems, [μ], ϕi=s)