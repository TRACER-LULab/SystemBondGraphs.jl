using BondGraphs
using ModelingToolkit
using DifferentialEquations
using MetaGraphs
using Graphs
## 
affinity_ATP_hydrolysis = 50000 / 8.314 / 310
S = Symbol
## Swap new node for old node and delete old node. 
function swap(BG, old, new)
    old = S(old)
    new = S(new)
    in_nodes = inneighbors(BG.graph, BG.graph[old, :name])
    out_nodes = outneighbors(BG.graph, BG.graph[old, :name])
    if !isempty(in_nodes)
        for i ∈ in_nodes
            add_edge!(BG.graph, i, BG.graph[new, :name])
        end
    end
    if !isempty(out_nodes)
        for i ∈ out_nodes
            add_edge!(BG.graph, BG.graph[new, :name], i)
        end
    end
    rem_vertex!(BG.graph, BG.graph[old, :name])
end
## Apply defaults from "old" to "new" element of the same type
function swap_defaults(BG, name, old, new)
    for key in collect(keys(BG[S(name*old)].defaults))
        BG[S(name*new)].defaults[key] = BG[S(name*old)].defaults[key]
    end
end
## Swap element for a bond element
function swap_bond(BG, elem)
    old_name = S(elem)
    new_name = S(elem*"_bond")
    add_Bond!(BG, new_name)
    swap(BG, elem, string(new_name))
    set_prop!(BG.graph, BG.graph[new_name, :name], :name, old_name)
end
## Kinase Factory
function kinase_factory!(BG, name; KM=1.0, KE=1.0)
    add_Se!(BG, S(name*"ATP"))
    add_Se!(BG, S(name*"ADP"))
    add_Ce!(BG, S(name*"E"))
    add_Ce!(BG, S(name*"C"))
    add_Ce!(BG, S(name*"MP"))
    add_Ce!(BG, S(name*"M"))
    for i ∈ 1:6
        add_Bond!(BG, S(name*"B$(i)"))
    end
    add_1J!(BG, Dict(
        S(name*"B1") => true,
        S(name*"B2") => false,
        S(name*"ATP") => true,
        S(name*"M") => true
        ), S(name*"J_M"))
    add_1J!(BG, Dict(
        S(name*"B5") => true,
        S(name*"MP") => false,
        S(name*"ADP") => false,
        S(name*"B6") => false
        ), S(name*"J_MP"))
    add_0J!(BG, Dict(
        S(name*"B3") => true, 
        S(name*"B4") => false,
        S(name*"C") => false
        ), S(name*"J_C"))
    add_0J!(BG, Dict(
        S(name*"B6") => true,
        S(name*"E") => false,
        S(name*"B1") => false
        ), S(name*"J_E"))
    add_Re!(BG, S(name*"B2"), S(name*"B3"), S(name*"Re1"))
    add_Re!(BG, S(name*"B4"), S(name*"B5"), S(name*"Re2"))
    kinase_parameters!(BG, name, KM=KM, KE=KE)
end
## Set Kinase Parameters
function kinase_parameters!(BG, name; KM=1.0, KE=1.0)
    ATP_potential = affinity_ATP_hydrolysis
    ADP_potential = 0.0
    a = 1000.0
    d = 150.0
    k = 150.0
    KC = exp(ATP_potential)*d/a*KM*KE
    r1 = d/KC
    r2 = k/KC
    BG[S(name*"E")].k = KE
    BG[S(name*"C")].k = KC
    BG[S(name*"M")].k = KM
    BG[S(name*"MP")].k = 0.0
    BG[S(name*"ATP")].Se = ATP_potential
    BG[S(name*"ADP")].Se = ADP_potential
    BG[S(name*"Re1")].r = r1
    BG[S(name*"Re2")].r = r2
end
## Phosphatase Factory
function phosphatase_factory!(BG, name; KMP=1.0, KE=1.0)
    add_Se!(BG, S(name*"P"))
    add_Ce!(BG, S(name*"E"))
    add_Ce!(BG, S(name*"MP"))
    add_Ce!(BG, S(name*"M"))
    add_Ce!(BG, S(name*"C"))
    for i ∈ 1:6
        add_Bond!(BG, S(name*"B$(i)"))
    end
    add_1J!(BG, Dict(
        S(name*"B1") => true,
        S(name*"B2") => false,
        S(name*"MP") => true
        ), S(name*"J_MP"))
    add_1J!(BG, Dict(
        S(name*"B5") => true,
        S(name*"B6") => false,
        S(name*"M") => false,
        S(name*"P") => false
        ), S(name*"J_M"))
    add_0J!(BG, Dict(
        S(name*"B3") => true,
        S(name*"C") => false,
        S(name*"B4") => false
        ), S(name*"J_C"))
    add_0J!(BG, Dict(
        S(name*"B6") => true,
        S(name*"B1") => false,
        S(name*"E") => false
        ), S(name*"J_E"))
    add_Re!(BG, S(name*"B2"), S(name*"B3"), S(name*"Re1"))
    add_Re!(BG, S(name*"B4"), S(name*"B5"), S(name*"Re2"))
    phosphatase_parameters!(BG, name, KMP=KMP, KE=KE)
end
## Phosphatase Parameters
function phosphatase_parameters!(BG, name; KMP=1.0, KE=1.0)
    P_potential = 0.0
    a = 1000.0
    d = 150.0
    k = 150.0
    KC = d/a*KMP*KE
    r1 = d/KC
    r2 = k/KC
    BG[S(name*"E")].k = KE
    BG[S(name*"C")].k = KC
    BG[S(name*"M")].k = 0.0
    BG[S(name*"MP")].k = KMP
    BG[S(name*"P")].Se = P_potential
    BG[S(name*"Re1")].r = r1
    BG[S(name*"Re2")].r = r2
end
## Phosphorylation Cycle
function phosphorylation_cycle_factory!(BG, name;KM=1.0, KMP=1.0, KKin=1.0)
    # Ports
    add_Ce!(BG, S(name*"M"))
    add_Ce!(BG, S(name*"MP"))
    add_Se!(BG, S(name*"ATP"))
    add_Se!(BG, S(name*"ADP"))
    add_Se!(BG, S(name*"P"))
    add_Ce!(BG, S(name*"Kin"))
    add_Ce!(BG, S(name*"Pho"))
    # Kinase and Phosphatase Reactions
    kinase_factory!(BG, name*"K_", KM=KM, KE=KKin)
    phosphatase_factory!(BG, name*"P_", KMP=KMP)
    phosphorylation_cycle_parameters!(BG, name, KMP=KMP, KM=KM, KKin=KKin)
    # Input and Output ports for cycle
    add_0J!(BG, Dict(
        S(name*"M") => true,
        S(name*"K_M") => false,
        S(name*"P_M") => true
        ), S(name*"J_M"))
    add_0J!(BG, Dict(
        S(name*"MP") => false,
        S(name*"K_MP") => true,
        S(name*"P_MP") => false
        ), S(name*"J_MP"))
    # Swap Species inputs for bond elements
    swap_bond(BG, name*"K_M")
    swap_bond(BG, name*"K_MP")
    swap_bond(BG, name*"P_MP")
    swap_bond(BG, name*"P_M")
    # Swap Components to create Cycle Ports
    swap(BG, name*"K_ATP", name*"ATP")
    swap(BG, name*"K_ADP", name*"ADP")
    swap(BG, name*"K_E", name*"Kin")
    swap(BG, name*"P_P", name*"P")
    swap(BG, name*"P_E", name*"Pho")
end
## Phosphorylation Cycle
function phosphorylation_cycle_parameters!(BG, name; KMP=1.0, KM=1.0, KKin=1.0)
    BG[S(name*"M")].k = KM
    BG[S(name*"MP")].k = KMP
    BG[S(name*"Kin")].k = KKin
    swap_defaults(BG, name, "K_ATP", "ATP")
    swap_defaults(BG, name, "K_ADP", "ADP")
    swap_defaults(BG, name, "K_E", "Kin")
    swap_defaults(BG, name, "P_P", "P")
    swap_defaults(BG, name, "P_E", "Pho") 
end
## MAPK Cascade Factory
function mapk_cascade_factory!(BG, name)
    species = [
        "MKKKK","MKKK","MKKKP","MKK","MKKP","MKKPP","MK","MKP","MKPP",
        "MAPKKKPase","MAPKKPase","MAPKPase"
        ]
    species_dict = Dict{String, Dict{Symbol, Bool}}()
    for s ∈ species
        add_Ce!(BG, S(name*s))
        species_dict[s] = Dict(S(name*s)=>true)
    end

    chemostat_dict = Dict{String, Dict{Symbol, Bool}}()
    chemostats = ["ATP", "ADP", "P"]
    for c ∈ chemostats
        add_Se!(BG, S(name*c))
        chemostat_dict[c] = Dict(S(name*c)=>true)
    end
    a = 1000.0
    d = 150.0
    k = 150.0
    D = (a*k/d)^2*exp(-affinity_ATP_hydrolysis)
    
    K_MKKKK = 1.0
    K_MKKK = 1.0
    K_MKK = 1.0
    K_MK = 1.0
    K_MKKKPase = 1.0
    K_MKKPase = 1.0
    K_MKPase = 1.0

    P_factor = sqrt(D)/k*d/a*exp(affinity_ATP_hydrolysis)
    K_MKKKP = K_MKKK*P_factor
    K_MKKP = K_MKK*P_factor
    K_MKKPP = K_MKKP*P_factor
    K_MKP = K_MK*P_factor
    K_MKPP = K_MKP*P_factor

    Kdict = Dict(
            "MKKKK"=> K_MKKKK,
            "MKKK"=> K_MKKK,
            "MKKKP"=> K_MKKKP,
            "MKK"=> K_MKK,
            "MKKP"=> K_MKKP,
            "MKKPP"=> K_MKKPP,
            "MK"=> K_MK,
            "MKP"=> K_MKP,
            "MKPP"=> K_MKPP,
            "MAPKKKPase"=> 1.0,
            "MAPKKPase"=> 1.0,
            "MAPKPase"=> 1.0,
            )
    phosphorylation_cycles = [
            Dict(:c => "cycle1_", :k => "MKKKK", :ph => "MAPKKKPase", :s => "MKKK", :p => "MKKKP"),
            Dict(:c => "cycle2_", :k => "MKKKP", :ph => "MAPKKPase",  :s => "MKK",  :p => "MKKP"),
            Dict(:c => "cycle3_", :k => "MKKKP", :ph => "MAPKKPase",  :s => "MKKP", :p => "MKKPP"),
            Dict(:c => "cycle4_", :k => "MKKPP", :ph => "MAPKPase",   :s => "MK",   :p => "MKP"),
            Dict(:c => "cycle5_", :k => "MKKPP", :ph => "MAPKPase",   :s => "MKP",  :p => "MKPP")
            ]
    for c ∈ phosphorylation_cycles
        phosphorylation_cycle_factory!(BG, c[:c]; KM=Kdict[c[:s]], KMP=Kdict[c[:p]], KKin=Kdict[c[:k]])
        if c[:c] == "cycle1_"
            for chemostat ∈ chemostats
                swap_defaults(BG, name, c[:c]*chemostat, chemostat)
            end
        end
        swap_bond(BG, name*c[:c]*"Kin")
        swap_bond(BG, name*c[:c]*"Pho")
        swap_bond(BG, name*c[:c]*"M")
        swap_bond(BG, name*c[:c]*"MP")
        swap_bond(BG, name*c[:c]*"ATP")
        swap_bond(BG, name*c[:c]*"ADP")
        swap_bond(BG, name*c[:c]*"P")
        species_dict[c[:k]][S(name*c[:c]*"Kin")]   = !true
        species_dict[c[:ph]][S(name*c[:c]*"Pho")]  = !true
        species_dict[c[:s]][S(name*c[:c]*"M")]     = !true
        species_dict[c[:p]][S(name*c[:c]*"MP")]    = !true
        chemostat_dict["ATP"][S(name*c[:c]*"ATP")] = !true
        chemostat_dict["ADP"][S(name*c[:c]*"ADP")] = !true
        chemostat_dict["P"][S(name*c[:c]*"P")]     = !true
    end
    for s ∈ species
        BG[S(name*s)].k = Kdict[s]
        add_0J!(BG, species_dict[s], S(name*"J_"*s))
    end
    for c ∈ chemostats 
        add_0J!(BG, chemostat_dict[c], S(name*"J_"*c))
    end
end
##
@parameters t
test = BondGraph(t) 
phosphorylation_cycle_factory!(test, "")
model = generate_model(test)
##
test = BondGraph(t)
mapk_cascade_factory!(test, "")
model = generate_model(test)
model = structural_simplify(model)
##
u0 = [
    test[:MKKKK].q      => 3e-5,
    test[:MKKK].q       => 3e-3,
    test[:MKK].q        => 1.2,
    test[:MK].q         => 1.2,
    test[:MAPKKKPase].q => 3e-4,
    test[:MAPKKPase].q  => 3e-4,
    test[:MAPKPase].q   => 1.2e-1
]
tspan = (0.0, 0.1)
prob = ODEProblem(model, u0, tspan, [])
sol = solve(prob)
## Plotting
import Pkg; Pkg.activate()
nodenames = map(i->string(test.graph.vprops[i][:name]), 1:nv(test.graph))
layout=(args...)->spring_layout(args...; C=70)
gplot(test.graph, nodelabel = nodenames, layout = layout, linetype="curve")
## 
using Symbolics, SymbolicUtils
##
explog = @rule exp(log(~x))=>~x
loglog = @rule 
@variables x_0, x_1, x_2, x_12
@rule(log(~x)+log(~y) => log(~x*~y))(log(x_2) + log(x_0))
##
file = open("equations.tex", "w") 
eqns = equations(model)
for i ∈ eachindex(eqns)
    write(file, latexify(eqns[i]))
end
close(file)