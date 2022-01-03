using BondGraphs
using ModelingToolkit
using DifferentialEquations
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
            S(name * "ATP") => false,
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
test = BioBondGraph(t, R = 1.0, T = 1.0)
mapk_cascade_factory!(test, "")
model = generate_model(test)
## Figure 5a
u0 = [
    test[:MKKKK].q => 3e-5,
    test[:MKKK].q => 3e-3,
    test[:MKK].q => 1.2,
    test[:MK].q => 1.2,
    test[:MAPKKKPase].q => 3e-4,
    test[:MAPKKPase].q => 3e-4,
    test[:MAPKPase].q => 1.2e-1
] |> Dict

tspan = (0.0, 100.0)
prob = ODEProblem(model, u0, tspan, [])
sol = solve(prob)
mkkkq = @. 100 * sol[test[:MKKKP].q] / sol[test[:MKKK].q][begin]
mkkq = @. 100 * sol[test[:MKKPP].q] / sol[test[:MKK].q][begin]
mkq = @. 100 * sol[test[:MKPP].q] / sol[test[:MK].q][begin]
plot(sol.t, [mkkkq, mkkq, mkq], ylims = (0, 100), ylabel = "Activation [%]", xlabel = "Time [s]", labels = ["MKKK" "MKK" "MK"], legend = :topleft, title = "Figure 5a")
## Figure 5b
init_mkkkk = 10 .^ range(-7, -1, 100)
res = similar(init_mkkkk, (3, length(init_mkkkk)))
u0 = Dict{Any,Float64}(states(model) .=> prob.u0)
for index in eachindex(init_mkkkk)
    u0[test[:MKKKK].q] = init_mkkkk[index]
    prob = remake(prob, u0 = ModelingToolkit.varmap_to_vars(u0, states(model)), tspan = (0.0, 1000.0))
    sol = solve(prob)
    res[:, index] .= (sol[test[:MKKKP].q][end], sol[test[:MKKPP].q][end], sol[test[:MKPP].q][end])
end
res[1, :] .= res[1, :] ./ maximum(res[1, :])
res[2, :] .= res[2, :] ./ maximum(res[2, :])
res[3, :] .= res[3, :] ./ maximum(res[3, :])
plot(init_mkkkk, res', xscale = :log, labels = ["MKPP" "MKKPP" "MKKKP"], legend = :topright, xlabel = "Input (MKKKK) Concentration (μM)", ylabel = "Normalized Concentration", title = "Figure 5b")
## Figure 5c
model.defaults[test[:ATP].Se] = model.defaults[test[:ATP].Se] * 0.8
p_new = ModelingToolkit.varmap_to_vars(model.defaults, parameters(model))
init_mkkkk = 10 .^ range(-7, 1, 100)
res = similar(init_mkkkk, (3, length(init_mkkkk)))
u0 = Dict{Any,Float64}(states(model) .=> prob.u0)
for index in eachindex(init_mkkkk)
    u0[test[:MKKKK].q] = init_mkkkk[index]
    prob = remake(prob, u0 = ModelingToolkit.varmap_to_vars(u0, states(model)), tspan = (0.0, 1000.0), p = p_new)
    sol = solve(prob)
    res[:, index] .= (sol[test[:MKKKP].q][end], sol[test[:MKKPP].q][end], sol[test[:MKPP].q][end])
end
res[1, :] .= res[1, :] ./ maximum(res[1, :])
res[2, :] .= res[2, :] ./ maximum(res[2, :])
res[3, :] .= res[3, :] ./ maximum(res[3, :])
plot(init_mkkkk, res', xscale = :log, labels = ["MKPP" "MKKPP" "MKKKP"], legend = :topleft, xlabel = "Input (MKKKK) Concentration (μM)", ylabel = "Normalized Concentration", title = "Figure 5c", xticks = 10.0.^(-4:2:0))
 