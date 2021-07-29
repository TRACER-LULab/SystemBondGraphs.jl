using DrWatson
@quickactivate "BondGraphModeling"
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##
using BondGraphs
using ModelingToolkit
using DifferentialEquations
## Setup Empty BondGraph
@variables t
visco = BondGraph(t);
## Inputs
add_Se!(visco, :σ₁)
add_Se!(visco, :σ₂)
add_Se!(visco, :εE)
## Create Connecting Bonds
add_Bond!(visco, :b2)
add_Bond!(visco, :b11)
add_Bond!(visco, :b5)
add_Bond!(visco, :b6)
add_Bond!(visco, :b8)
add_Bond!(visco, :b9)
## Multiport C Elements
# Stress-relationship
function ϕi(𝐞, 𝐪, params)
    λ₁, λ₂ = 𝐪
    λ₃ = 1 / λ₁ / λ₂
    μ = params
    σ1 = μ * (λ₁^2 - λ₃^2)
    σ2 = μ * (λ₂^2 - λ₃^2)
    return [σ1; σ2]
end 
@parameters μ
# elements
elems = [:b5 => false, :b8 => false]
add_C_multiport!(visco, elems, [μ], :Cα, ϕi = ϕi)
elems = [:b6 => false, :b9 => false]
add_C_multiport!(visco, elems, [μ], :Cβ, ϕi = ϕi)
## Add Dampers
add_R!(visco, :R1)
add_R!(visco, :R2)
## Add 1-junctions
add_1J!(visco, Dict([
    :σ₁ => false, 
    :εE => false,
    :b5 => true,
    :b2 => true
    ]), :J1_1)
add_1J!(visco, Dict([
    :σ₂ => false, 
    :εE => false,
    :b8 => true,
    :b11 => true
    ]), :J1_2)
## Add 0-Junctions
add_0J!(visco, Dict([
    :b2 => false,
    :R1 => true,
    :b6 => true
    ]), :J0_1)
add_0J!(visco, Dict([
    :b11 => false,
    :R2 => true,
    :b9 => true
    ]), :J0_2)
## Generate the model
generate_model!(visco)

