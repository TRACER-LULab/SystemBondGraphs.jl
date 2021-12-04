using BondGraphs
using ModelingToolkit
using DifferentialEquations
using MetaGraphs
using Graphs
##
affinity_ATP_hydrolysis = 50000 / 8.314 / 310
S = Symbol 
##
include("helpers.jl")
include("kinase_factory.jl")
include("phosphatase_factory.jl")
include("phosphorylation_cycle.jl")
## MAPK Cascade Factory
function mapk_cascade_factory!(BG, name)
    species = [
        "MKKKK", "MKKK", "MKKKP", "MKK", "MKKP", "MKKPP", "MK", "MKP", "MKPP",
        "MAPKKKPase", "MAPKKPase", "MAPKPase"
    ]
    species_dict = Dict{String,Dict{Symbol,Bool}}()
    for s ∈ species
        add_Ce!(BG, S(name * s))
        species_dict[s] = Dict(S(name * s) => !true)
    end

    chemostat_dict = Dict{String,Dict{Symbol,Bool}}()
    chemostats = ["ATP", "ADP", "P"]
    for c ∈ chemostats
        add_Se!(BG, S(name * c))
        chemostat_dict[c] = Dict(S(name * c) => true)
    end
    a = 1000.0
    d = 150.0
    k = 150.0
    D = (a * k / d)^2 * exp(-affinity_ATP_hydrolysis)

    K_MKKKK = 1.0
    K_MKKK = 1.0
    K_MKK = 1.0
    K_MK = 1.0
    K_MKKKPase = 1.0
    K_MKKPase = 1.0
    K_MKPase = 1.0

    P_factor = sqrt(D) / k * d / a * exp(affinity_ATP_hydrolysis)
    K_MKKKP = K_MKKK * P_factor
    K_MKKP = K_MKK * P_factor
    K_MKKPP = K_MKKP * P_factor
    K_MKP = K_MK * P_factor
    K_MKPP = K_MKP * P_factor

    Kdict = Dict(
        "MKKKK" => K_MKKKK,
        "MKKK" => K_MKKK,
        "MKKKP" => K_MKKKP,
        "MKK" => K_MKK,
        "MKKP" => K_MKKP,
        "MKKPP" => K_MKKPP,
        "MK" => K_MK,
        "MKP" => K_MKP,
        "MKPP" => K_MKPP,
        "MAPKKKPase" => 1.0,
        "MAPKKPase" => 1.0,
        "MAPKPase" => 1.0,
    )
    phosphorylation_cycles = [
        Dict(:c => "cycle1_", :k => "MKKKK", :ph => "MAPKKKPase", :s => "MKKK", :p => "MKKKP"),
        Dict(:c => "cycle2_", :k => "MKKKP", :ph => "MAPKKPase", :s => "MKK", :p => "MKKP"),
        Dict(:c => "cycle3_", :k => "MKKKP", :ph => "MAPKKPase", :s => "MKKP", :p => "MKKPP"),
        Dict(:c => "cycle4_", :k => "MKKPP", :ph => "MAPKPase", :s => "MK", :p => "MKP"),
        Dict(:c => "cycle5_", :k => "MKKPP", :ph => "MAPKPase", :s => "MKP", :p => "MKPP")
    ]
    for c ∈ phosphorylation_cycles
        phosphorylation_cycle_factory!(BG, c[:c]; KM = Kdict[c[:s]], KMP = Kdict[c[:p]], KKin = Kdict[c[:k]])
        if c[:c] == "cycle1_"
            for chemostat ∈ chemostats
                swap_defaults(BG, name, c[:c] * chemostat, chemostat)
            end
        end
        swap_bond(BG, name * c[:c] * "Kin")
        swap_bond(BG, name * c[:c] * "Pho")
        swap_bond(BG, name * c[:c] * "M")
        swap_bond(BG, name * c[:c] * "MP")
        swap_bond(BG, name * c[:c] * "ATP")
        swap_bond(BG, name * c[:c] * "ADP")
        swap_bond(BG, name * c[:c] * "P")
    end
    # Connect MKKKK
    add_0J!(BG, Dict(
            S(name * "MKKKK") => false,
            S(name * "cycle1_Kin") => true
        ), S(name * "J_MKKKK"))
    # Connect MKKK
    add_0J!(BG, Dict(
            S(name * "MKKK") => false,
            S(name * "cycle1_M") => true,
        ), S(name * "J_MKKK"))
    # Connect MKKKPase
    add_0J!(BG, Dict(
            S(name * "MAPKKKPase") => false,
            S(name * "cycle1_Pho") => true
        ), S(name * "J_MAPKKKPase"))
    # Connect MKKKP
    add_0J!(BG, Dict(
            S(name * "MKKKP") => false,
            S(name * "cycle1_MP") => true,
            S(name * "cycle2_Kin") => true,
            S(name * "cycle3_Kin") => true
        ), S(name * "J_MKKKP"))
    # Connect MKK
    add_0J!(BG, Dict(
            S(name * "MKK") => false,
            S(name * "cycle2_M") => true
        ), S(name * "J_MKK"))
    # Connect MKKP
    add_0J!(BG, Dict(
            S(name * "MKKP") => false,
            S(name * "cycle2_MP") => true,
            S(name * "cycle3_M") => true
        ), S(name * "J_MKKP"))
    # Connect MAPKK-Pase
    add_0J!(BG, Dict(
            S(name * "MAPKKPase") => false,
            S(name * "cycle2_Pho") => true,
            S(name * "cycle3_Pho") => true
        ), S(name * "J_MAPKKPase"))
    # Connect MAPKKPP
    add_0J!(BG, Dict(
            S(name * "MKKPP") => false,
            S(name * "cycle3_MP") => true,
            S(name * "cycle4_Kin") => true,
            S(name * "cycle5_Kin") => true
        ), S(name * "J_MKKPP"))
    # Connect MAPK
    add_0J!(BG, Dict(
            S(name * "MK") => false,
            S(name * "cycle4_M") => true
        ), S(name * "J_MK"))
    # Connect MKP
    add_0J!(BG, Dict(
            S(name * "MKP") => false,
            S(name * "cycle4_MP") => true,
            S(name * "cycle5_M") => true
        ), S(name * "J_MKP"))
    # Connect MAPKPase    
    add_0J!(BG, Dict(
            S(name * "MAPKPase") => false,
            S(name * "cycle4_Pho") => true,
            S(name * "cycle5_Pho") => true,
        ), S(name * "J_MAPKPase"))
    # Connect MKPP
    add_0J!(BG, Dict(
            S(name * "MKPP") => false,
            S(name * "cycle5_MP") => true
        ), S(name * "J_MKPP"))
    # Connect ATP
    add_0J!(BG, Dict(
            S(name * "ATP") => true,
            S(name * "cycle1_ATP") => false,
            S(name * "cycle2_ATP") => false,
            S(name * "cycle3_ATP") => false,
            S(name * "cycle4_ATP") => false,
            S(name * "cycle5_ATP") => false,
        ), S(name * "J_ATP"))
    # Connect ADP
    add_0J!(BG, Dict(
            S(name * "ADP") => false,
            S(name * "cycle1_ADP") => true,
            S(name * "cycle2_ADP") => true,
            S(name * "cycle3_ADP") => true,
            S(name * "cycle4_ADP") => true,
            S(name * "cycle5_ADP") => true,
        ), S(name * "J_ADP"))
    # Connect P
    add_0J!(BG, Dict(
            S(name * "P") => false,
            S(name * "cycle1_P") => true,
            S(name * "cycle2_P") => true,
            S(name * "cycle3_P") => true,
            S(name * "cycle4_P") => true,
            S(name * "cycle5_P") => true,
        ), S(name * "J_P"))
    for s ∈ species
        BG[S(name * s)].k = Kdict[s]
    end
end
##
@parameters t
test = BioBondGraph(t, R=1.0, T=1.0)
mapk_cascade_factory!(test, "")
model = generate_model(test)
model = structural_simplify(model, simplify = true)
##
using SymbolicUtils
RW = SymbolicUtils.Rewriters
r1 = @acrule log(~x) + log(~y) => log((~x) * (~y))
r2 = @rule log(~x) - log(~y) => log((~x) / (~y))
r3 = @rule (~x) * log(~y) => log((~y)^(~x))
r4 = @rule exp(log(~x)) => ~x
r5 = @acrule exp((~x) + (~y)) => exp(~x) * exp(~y)
rw1 = RW.Fixpoint(RW.Chain([r1, r2, r3, r4, r5]))
rw2 = RW.Prewalk(RW.Chain([r1, r2, r3, r4, r5]))
rw3 = RW.Postwalk(RW.Chain([r1, r2, r3, r4, r5]))
eqns = equations(model)
defaults = model.defaults
params = parameters(model)
default_params = Dict(map(x -> x => defaults[x], params))
eqns = map(eqn -> substitute(eqn, default_params) |> expand |> simplify, eqns)
for i ∈ eachindex(eqns)
    eqns[i] = eqns[i].lhs ~ eqns[i].rhs |> rw3 |> rw2 |> rw1 |> expand
end
eqns[1]
##
sys = ODESystem(eqns, name = :model, defaults = defaults)
model = structural_simplify(sys, simplify = true)
##
## ADDED RULES
# @rule(sum(log, ~~xs) => log(*(~~xs)))
# @rule(exp(log(~x)) => ~x)
# @rule(~x*log(~y) => log((~y) ^ (~x)))
# @rule(exp(+(~~xs))=>prod(exp, ~~xs))
u0 = [
    test[:MKKKK].q => 3e-5,
    test[:MKKK].q => 3e-3,
    test[:MKK].q => 1.2,
    test[:MK].q => 1.2,
    test[:MAPKKKPase].q => 3e-4,
    test[:MAPKKPase].q => 3e-4,
    test[:MAPKPase].q => 1.2e-1
]
tspan = (0.0, 100.0)
prob = ODEProblem(model, u0, tspan, [])
sol = solve(prob)
## Plotting
import Pkg;
Pkg.activate();
nodenames = map(i -> string(test.graph.vprops[i][:name]), 1:nv(test.graph))
layout = (args...) -> spring_layout(args...; C = 70)
gplot(test.graph, nodelabel = nodenames, layout = layout, linetype = "curve")
## 
using Symbolics, SymbolicUtils
##
explog = @rule exp(log(~x)) => ~x
loglog = @rule
@syms x y z
acr1 = @acrule(log(~x) + log(~y) => log(~x * ~y))
acr1(log(x) + log(y))
##
@parameters t
@variables e_0(t) e_1(t) e_2(t) e_3(t) e_4(t) e_5(t) e_6(t)
@variables f_0(t) f_1(t) f_2(t) f_3(t) f_4(t) f_5(t) f_6(t)
@variables x_0(t) x_1(t) x_2(t) x_3(t) x_4(t) x_5(t) x_6(t)
D = Differential(t)
eqn = [D(x_0) ~ -300 * x_0 + 3.75624385638826e-6 * (exp(e_0 + e_2 + e_4) - exp(e_0 + e_3 + e_5)),
    D(x_1) ~ -300 * x_1 + 1000 * (exp(e_1 + e_3) - exp(e_1 + e_2 + e_6)),
    D(x_2) ~ -1 * (150 * x_0 + 150 * x_1 - 3.75624385638826e-6 * exp(e_0 + e_2 + e_4) - 1000 * exp(e_1 + e_2 + e_6)),
    D(x_3) ~ -1 * (150 * x_0 + 150 * x_1 - 1000 * exp(e_1 + e_3) - 3.75624385638826e-6 * exp(e_0 + e_3 + e_5)),
    D(x_4) ~ -1 * (150 * x_0 - 3.75624385638826e-6 * exp(e_0 + e_2 + e_4)),
    D(x_5) ~ -1 * (150 * x_0 - 3.75624385638826e-6 * exp(e_0 + e_3 + e_5)),
    D(x_6) ~ -1 * (150 * x_1 - 1000 * exp(e_1 + e_2 + e_6)),
    e_0 ~ log(x_0),
    e_1 ~ log(x_1),
    e_2 ~ log(x_2),
    e_3 ~ log(x_3),
    e_4 ~ log(x_4),
    e_5 ~ log(x_5),
    e_6 ~ log(x_6),
]
@named sys = ODESystem(eqn)
sys = structural_simplify(sys)
##
file = open("equations.tex", "w")
eqns = equations(model)
for i ∈ eachindex(eqns)
    write(file, latexify(eqns[i]))
end
close(file)


cycle1_K_Re1₊r * cycle1_K_C₊k * cycle1_K_C₊q(t)
+
cycle1_K_Re2₊r * MKKKK₊k * MKKKK₊q(t) * MKKKP₊k * MKKKP₊q(t)
-
cycle1_K_Re1₊r * MKKKK₊k * MKKKK₊q(t) * MKKK₊k * MKKK₊q(t) * exp(ATP₊Se(t))
-
cycle1_K_Re2₊r * cycle1_K_C₊k * cycle1_K_C₊q(t)