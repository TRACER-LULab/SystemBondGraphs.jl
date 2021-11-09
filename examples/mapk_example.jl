using BondGraphs
using ModelingToolkit
using DifferentialEquations
using MetaGraphs
using Graphs
## 
affinity_ATP_hydrolysis = 50000 / 8.314 / 310

function swap(BG, input, output)
    input = Symbol(input)
    output = Symbol(output)
    in_nodes = inneighbors(BG.graph, BG.graph[input, :name])
    out_nodes = outneighbors(BG.graph, BG.graph[input, :name])
    if !isempty(in_nodes)
        for i ∈ in_nodes
            add_edge!(BG.graph, i, BG.graph[output, :name])
        end
    end
    if !isempty(out_nodes)
        for i ∈ out_nodes
            add_edge!(BG.graph, BG.graph[output, :name], i)
        end
    end
    rem_vertex!(BG.graph, BG.graph[input, :name])
end

function connector(BG, input, output)
    set_prop!(BG.graph, BG.graph[Symbol(input), :name], :name, Symbol(output))
end

function kinase_factory!(BG::BondGraph, name::String; K_M = 1.0, K_E = 1.0)
#
    add_Bond!(BG, Symbol(name * "ATP"))
    add_Bond!(BG, Symbol(name * "ADP"))
    add_Bond!(BG, Symbol(name * "E"))
    add_Bond!(BG, Symbol(name * "M"))
    add_Bond!(BG, Symbol(name * "MP"))
    add_Bond!(BG, Symbol(name * "C"))

    for i ∈ 1:6
        add_Bond!(BG, Symbol(name * "B$i"))
    end

    add_1J!(BG, Dict(
        Symbol(name * "M") => true, 
        Symbol(name * "ATP") => true,
        Symbol(name * "B1") => true,
        Symbol(name * "B2") => false
        ),
        Symbol(name * "J11"))

    add_Re!(BG, Symbol(name * "B2"), Symbol(name * "B3"), Symbol(name * "re1"))

    add_0J!(BG, Dict(
        Symbol(name * "B3") => true, 
        Symbol(name * "B4") => false,
        Symbol(name * "C") => false
        ), Symbol(name * "J0C"))

    add_Re!(BG, Symbol(name * "B4"), Symbol(name * "B5"), Symbol(name * "re2"))

    add_1J!(BG, Dict(
        Symbol(name * "B5") => true,
        Symbol(name * "MP") => false,
        Symbol(name * "ADP") => false,
        Symbol(name * "B6") => false
        ), Symbol(name * "J12"))

    add_0J!(BG, Dict(
        Symbol(name * "E") => false,
        Symbol(name * "B6") => true,
        Symbol(name * "B1") => false
        ), Symbol(name * "J0E"))
#
    ATP_potential = affinity_ATP_hydrolysis
    ADP_potential = 0.0
    a = 1000.0
    d = 150.0
    k = 150.0
    K_C = exp(ATP_potential) * d / a * K_M * K_E
    r1 = d / K_C
    r2 = k / K_C
    # BG[Symbol(name * "E")].k = K_E
    # BG[Symbol(name * "C")].k = K_C
    # BG[Symbol(name * "M")].k = K_M
    # BG[Symbol(name * "MP")].k = 0.0
    # BG[Symbol(name * "ATP")].Se = ATP_potential
    # BG[Symbol(name * "ADP")].Se = ADP_potential
    # BG[Symbol(name * "re1")].r = r1
    # BG[Symbol(name * "re2")].r = r2
end
## Test Kinase Factory
    # kf = BondGraph(t)
    # kinase_factory!(kf, "")
    # model = generate_model(kf)
    # structural_simplify(model)
##
function phosphatase_factory!(BG::BondGraph, name::String; K_MP = 1.0, K_E = 1.0)
#    
    add_Bond!(BG, Symbol(name * "MP"))
    add_Bond!(BG, Symbol(name * "E"))
    add_Bond!(BG, Symbol(name * "M"))
    add_Bond!(BG, Symbol(name * "C"))
    add_Bond!(BG, Symbol(name * "Pi"))

    for i ∈ 1:6
        add_Bond!(BG, Symbol(name * "B$i"))
    end

    add_1J!(BG, Dict(
        Symbol(name * "MP") => true,     
        Symbol(name * "B1") => true,
        Symbol(name * "B2") => false
        ),
        Symbol(name * "J11"))

    add_Re!(BG, Symbol(name * "B2"), Symbol(name * "B3"), Symbol(name * "re1"))

    add_0J!(BG, Dict(
        Symbol(name * "B3") => true, 
        Symbol(name * "B4") => false,
        Symbol(name * "C") => false
        ), Symbol(name * "J0C"))

    add_Re!(BG, Symbol(name * "B4"), Symbol(name * "B5"), Symbol(name * "re2"))

    add_1J!(BG, Dict(
        Symbol(name * "B5") => true,
        Symbol(name * "M") => false,
        Symbol(name * "Pi") => false,
        Symbol(name * "B6") => false
        ), Symbol(name * "J12"))

    add_0J!(BG, Dict(
        Symbol(name * "E") => false,
        Symbol(name * "B6") => true,
        Symbol(name * "B1") => false
        ), Symbol(name * "J0E"))
#
    P_potential = 0.0
    a = 1000.0
    d = 150.0
    k = 150.0
    K_C = d / a * K_MP * K_E
    r1 = d / K_C
    r2 = k / K_C
    # BG[Symbol(name * "E")].k = K_E
    # BG[Symbol(name * "C")].k = K_C
    # BG[Symbol(name * "M")].k = 0.0
    # BG[Symbol(name * "MP")].k = K_MP
    # BG[Symbol(name * "Pi")].Se = P_potential
    # BG[Symbol(name * "re1")].r = r1
    # BG[Symbol(name * "re2")].r = r2
end
## Test Phosphatase Factor
    # pf = BondGraph(t)
    # phosphatase_factory!(pf, "")
    # model = generate_model(pf)
    # structural_simplify(model)
##
function phosphorylation_cycle_factory!(BG::BondGraph, name::String; K_M = 1.0, K_MP = 1.0, K_Kin = 1.0)
#    Create K and P factories
    kinase_factory!(BG, name * "K_", K_M = K_M, K_E = K_Kin)
    phosphatase_factory!(BG, name * "P_", K_MP = K_MP)
#    Create Species and Chemostat inputs
    add_Bond!(BG, Symbol(name * "M"))
    add_Bond!(BG, Symbol(name * "MP"))
    # add_Bond!(BG, Symbol(name * "ATP"))
    # add_Bond!(BG, Symbol(name * "ADP"))
    # add_Bond!(BG, Symbol(name * "Pi"))
    # add_Bond!(BG, Symbol(name * "Kin"))
    # add_Bond!(BG, Symbol(name * "Pho"))
#    Connect Everything
    add_0J!(BG, Dict(
        Symbol(name * "M") => true,
        Symbol(name * "K_M") => false,
        Symbol(name * "P_M") => true
        ),
        Symbol(name * "J0M"))

    add_0J!(BG, Dict(
        Symbol(name * "MP") => false,
        Symbol(name * "K_MP") => true,
        Symbol(name * "P_MP") => false
        ),
        Symbol(name * "J0MP"))

    connector(BG, name * "K_ATP", name * "ATP")
    connector(BG, name * "K_E",   name * "Kin")
    connector(BG, name * "K_ADP", name * "ADP")
    connector(BG, name * "P_Pi",  name * "Pi")
    connector(BG, name * "P_E",   name * "Pho")
# Set Parameters
    # BG[Symbol(name * "M")].k = K_M
    # BG[Symbol(name * "MP")].k = K_MP
    # BG[Symbol(name * "Kin")].k = K_Kin
end
##
# pcf = BondGraph(t)
# phosphorylation_cycle_factory!(pcf, "")
# model = generate_model(pcf)
# structural_simplify(model)
##
function mapk_factory!(BG::BondGraph, name::String)
#
    a = 1000.0
    d = 150.0
    k = 150.0
    D = (a*k/d)^2*exp(-affinity_ATP_hydrolysis)
    K_M4K = 1.0
    K_M3K = 1.0
    K_M2K = 1.0
    K_MK = 1.0
    K_M3KPase = 1.0
    K_M2KPase = 1.0
    K_MKPase = 1.0

    P_factor = (sqrt(D)/k)*(d/a)*exp(affinity_ATP_hydrolysis)
    K_M3KP = K_M3K*P_factor
    K_M2KP = K_M2K*P_factor
    K_M2KPP = K_M2KP*P_factor
    K_MKP = K_MK*P_factor
    K_MKPP = K_MKP*P_factor

#
    Kdict = Dict(
        "M4K"=> K_M4K,
        "M3K"=> K_M3K,
        "M3KP"=> K_M3KP,
        "M2K"=> K_M2K,
        "M2KP"=> K_M2KP,
        "M2KPP"=> K_M2KPP,
        "MK"=> K_MK,
        "MKP"=> K_MKP,
        "MKPP"=> K_MKPP,
        "MAP3K-Pase"=> 1.0,
        "MAP2K-Pase"=> 1.0,
        "MAPK-Pase"=> 1.0,
    )
#
    add_Ce!(BG, Symbol(name * "MAP4K"))
    add_Ce!(BG, Symbol(name * "MAP3K"))
    add_Ce!(BG, Symbol(name * "MAP3KPase"))
    add_Ce!(BG, Symbol(name * "MAP3KP"))
    add_Ce!(BG, Symbol(name * "MAP2K"))
    add_Ce!(BG, Symbol(name * "MAP2KPase"))
    add_Ce!(BG, Symbol(name * "MAP2KP"))
    add_Ce!(BG, Symbol(name * "MAP2KPP"))
    add_Ce!(BG, Symbol(name * "MAPK"))
    add_Ce!(BG, Symbol(name * "MAPKPase"))
    add_Ce!(BG, Symbol(name * "MAPKP"))
    add_Ce!(BG, Symbol(name * "MAPKPP"))

# 
    Phosphorylation_cycles = [
        Dict(:name=>"cycle1_", :kin=>"M4K",   :pho=>"MAP3KPase", :s=>"M3K",  :p=>"M3KP"),
        Dict(:name=>"cycle2_", :kin=>"M3KP",  :pho=>"MAP2KPase", :s=>"M2K",  :p=>"M2KP"),
        Dict(:name=>"cycle3_", :kin=>"M3KP",  :pho=>"MAP2KPase", :s=>"M2KP", :p=>"M2KPP"),
        Dict(:name=>"cycle4_", :kin=>"M2KPP", :pho=>"MAPKPase",  :s=>"MK",   :p=>"MKP"),
        Dict(:name=>"cycle5_", :kin=>"M2KPP", :pho=>"MAPKPase",  :s=>"MKP",  :p=>"MKPP"),
    ]
    # Create Phosphorylation Cycles
    for cycle in Phosphorylation_cycles
        phosphorylation_cycle_factory!(BG, name * cycle[:name], K_M = Kdict[cycle[:s]], K_MP = Kdict[cycle[:p]], K_Kin=Kdict[cycle[:kin]])
    end

    # Direct Connections from Species to Cyles
    swap(BG, name * "cycle1_M",   name * "MAP3K")
    swap(BG, name * "cycle1_Pho", name * "MAP3KPase")
    swap(BG, name * "cycle1_Kin", name * "MAP4K")
    swap(BG, name * "cycle2_M",   name * "MAP2K")
    swap(BG, name * "cycle4_M",   name * "MAPK")
    swap(BG, name * "cycle5_MP",  name * "MAPKPP")
    
    # Create 0-Junction Connections
    add_0J!(BG, Dict(
        Symbol(name * "cycle1_MP") => true,
        Symbol(name * "cycle3_Kin") => true,
        ), Symbol(name * "JMAP2KP"))
    add_0J!(BG, Dict(
        Symbol(name * "cycle2_Pho") => false,
        Symbol(name * "cycle3_Pho") => false,
        Symbol(name * "MAP2KPase") => false
        ), Symbol(name * "JMAP2KPase"))
    add_0J!(BG, Dict(
        Symbol(name * "cycle3_MP") => true, 
        Symbol(name * "cycle4_Kin") => false,
        Symbol(name * "cycle5_Kin") => false,
        Symbol(name * "MAP2KPP") => false
        ), Symbol(name * "JMAP2KPP"))
    add_0J!(BG, Dict(
        Symbol(name * "cycle4_MP") => true,
        Symbol(name * "cycle5_M") => true,
        Symbol(name * "MAPKP") => false
        ), Symbol(name * "JMAPKP"))
    add_0J!(BG, Dict(
        Symbol(name * "cycle4_Pho") => false,
        Symbol(name * "cycle5_Pho") => false,
        Symbol(name * "MAPKPase") => false
        ), Symbol(name * "JMAPKPase"))
    # Add Species for ATP, Pi, ADP
    add_Se!(BG, Symbol(name * "ATP"))
    add_Se!(BG, Symbol(name * "ADP"))
    add_Se!(BG, Symbol(name * "Pi"))

    atp_dict = Dict(Symbol(name * "ATP") => false)
    adp_dict = Dict(Symbol(name * "ADP") => false)
    pi_dict = Dict(Symbol(name * "Pi") => false)

    for i ∈ 1:5
        atp_dict[Symbol(name * "cycle$(i)_ATP")] = true
        adp_dict[Symbol(name * "cycle$(i)_ADP")] = true
        pi_dict[Symbol(name * "cycle$(i)_Pi")] = true
    end

    add_0J!(BG, atp_dict, Symbol(name * "JATP"))
    add_0J!(BG, adp_dict, Symbol(name * "JADP"))
    add_0J!(BG, pi_dict, Symbol(name * "JPi"))
end

## Create MAPK Cascade

@parameters t
mapk = BondGraph(t)
mapk_factory!(mapk, "")
##
model = generate_model(mapk)
model = structural_simplify(model)

## 
p = [
    # Good
    mapk[:MAP4K].k        => 1.0,
    mapk[:MAP3K].k        => 1.0,
    mapk[:MAP3KP].k       => 16316.353672239,
    mapk[:MAP2K].k        => 1.0,
    mapk[:MAP2KP].k       => 16316.353672239,
    mapk[:MAP2KPP].k      => 266223397.1575871,
    mapk[:MAPK].k         => 1.0,
    mapk[:MAPKP].k        => 16316.353672239,
    mapk[:MAPKPP].k       => 2.6622e8,
    mapk[:MAP3KPase].k    => 1.0,
    mapk[:MAP2KPase].k    => 1.0,
    mapk[:MAPKPase].k     => 1.0,
    # Good
    mapk[:cycle1_K_C].k   => 3.9934e7,
    mapk[:cycle1_P_C].k   => 2.4475e3,
    mapk[:cycle1_K_re1].r => 3.7562e-06,
    mapk[:cycle1_K_re2].r => 3.7562e-06,
    mapk[:cycle1_P_re1].r => 6.1288e-02,
    mapk[:cycle1_P_re2].r => 6.1288e-02,
    # Good
    mapk[:cycle2_K_C].k   => 6.5157e11,
    mapk[:cycle2_P_C].k   => 2.4475e3,
    mapk[:cycle2_K_re1].r => 2.3021e-10,
    mapk[:cycle2_K_re2].r => 2.3021e-10,
    mapk[:cycle2_P_re1].r => 6.1288e-02,
    mapk[:cycle2_P_re2].r => 6.1288e-02,
    # Good
    mapk[:cycle3_K_C].k   => 1.0631e16,
    mapk[:cycle3_P_C].k   => 3.9934e7,
    mapk[:cycle3_K_re1].r => 1.4109e-14,
    mapk[:cycle3_K_re2].r => 1.4109e-14,
    mapk[:cycle3_P_re1].r => 3.7562e-06,
    mapk[:cycle3_P_re2].r => 3.7562e-06,
    # Good
    mapk[:cycle4_K_C].k   => 1.0631e16,
    mapk[:cycle4_P_C].k   => 2.4475e3,
    mapk[:cycle4_K_re1].r => 1.4109e-14,
    mapk[:cycle4_K_re2].r => 1.4109e-14,
    mapk[:cycle4_P_re1].r => 6.1288e-02,
    mapk[:cycle4_P_re2].r => 6.1288E-02,
    # Good
    mapk[:cycle5_K_C].k   => 1.7346E+20,
    mapk[:cycle5_P_C].k   => 3.9934E+07,
    mapk[:cycle5_K_re1].r => 8.6474E-19,
    mapk[:cycle5_K_re2].r => 8.6474E-19,
    mapk[:cycle5_P_re1].r => 3.7562E-06,
    mapk[:cycle5_P_re2].r => 3.7562E-06,
    # Good
    mapk[:ATP].Se         => 19.4,
    mapk[:ADP].Se         => 0.0,
    mapk[:Pi].Se          => 0.0,
    ]
u0 = [
    mapk[:MAP4K].q => 0.03e-3,
    mapk[:MAP3K].q => 3e-3,
    mapk[:MAP2K].q => 1.2,
    mapk[:MAPK].q  => 1.2,
    mapk[:MAP3KPase].q => 3e-4,
    mapk[:MAP2KPase].q => 3e-4,
    mapk[:MAPKPase].q  => 0.12
    ]
prob = ODAEProblem(model, u0, (0.0, 100.0), p, jac = true, sparse = true)
sol = solve(prob, Rodas4())
# plot(sol, legend = false)