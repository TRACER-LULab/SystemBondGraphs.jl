using BondGraphs
using ModelingToolkit
using DifferentialEquations
using Plots
using Symbolics
using SymbolicUtils
# using Unitful
# using UnitfulRecipes
# using MonteCarloMeasurements
# using Latexify

@variables t
@parameters L x
rod = BondGraph(t)

n = 3
ρ = 7850.0
A = 1.0
L = 1.0
E = 210e9 # Elastic modulus of ASTM A228
Y(n, x) = sin((2 * n - 1) * π / 2 * x / L)
for i ∈ 1:n
    ωm = √(E / ρ) * (2 * i - 1) / L * π / 2
    mm = ρ * A * L / 2
    add_C!(rod, Symbol(string("C_", i)))
    rod.elements[Symbol(string("C_", i))].sys.C = mm * ωm^2
    add_I!(rod, Symbol(string("I_", i)))
    rod.elements[Symbol(string("I_", i))].sys.I = mm
    add_Bond!(rod, Symbol(string("B_", i)))
    add_1J!(rod, Dict([
        Symbol(string("C_", i)) => true,
        Symbol(string("I_", i)) => true,
        Symbol(string("B_", i)) => false       
        ]),
        Symbol(string("J1_", i)))
    add_Bond!(rod, Symbol(string("B_", i + n)))
    @parameters 
    add_TF!(rod, Sym{Num}(Symbol("m_" * string(i))), Dict([
        Symbol(string("B_", i)) => true,
        Symbol(string("B_", i + n)) => false
        ]),
        Symbol(string("TF_", i)))
end
add_Bond!(rod, :B1)
add_0J!(rod, Dict([map(i -> Symbol(string("B_", i + n)) => true, 1:n)...; :B1 => false]), :J0_F)
add_Bond!(rod, :B2)
add_1J!(rod, Dict([
    :B1 => true,
    :B2 => true
    ]),
    :J_nF)
add_Bond!(rod, :B3)
add_Bond!(rod, :B4)
add_0J!(rod, Dict([
    :B2 => false,
    :B3 => true,
    :B4 => true
    ]),
    :J0_1)
add_R!(rod, :R_b)
add_C!(rod, :C_k)
add_1J!(rod, Dict([
    :B3 => false,
    :R_b => true,
    :C_k => true
    ]),
    :J1_C)
add_I!(rod, :I_m)
F(t, p) = p[1] * sin(p[2] * t)
@register F(t, p)
@parameters α ω
add_Se!(rod, F, [α, ω], :S_e)
add_1J!(rod, Dict([
    :B4 => true,
    :I_m => true,
    :S_e => false
    ]),
    :J1_in)

generate_model!(rod)
simplify_model!(rod)

u0 = Dict(map(x -> x => 0.0, states(rod.model)))
u0 = ArrayPartition(map(x -> [u0[x]], states(rod.model))...)
p = [
    rod.elements[:I_m].sys.I => 1.0,
    rod.elements[:R_b].sys.R => 1.0,
    rod.elements[:C_k].sys.C => 0.1,
    rod.elements[:S_e].sys.α => 1.0,
    rod.elements[:S_e].sys.ω => 2.0,
] |> Dict
for i ∈ 1:n
    ωm = √(E / ρ) * (2 * i - 1) / L * π / 2
    mm = ρ * A * L / 2
    p[rod.elements[Symbol(string("C_", i))].sys.C] = mm * ωm^2
    p[rod.elements[Symbol(string("I_", i))].sys.I] = mm
end
ps = ArrayPartition(map(x -> [p[x]], parameters(rod.model))...)

tspan = (0.0, 10.0)

func = ODEFunction(rod.model)
prob = ODEProblem(func, u0, tspan, ps)
sol = solve(prob)
plot(sol)
