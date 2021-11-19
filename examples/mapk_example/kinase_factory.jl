## Kinase Factory
function kinase_factory!(BG, name; KM = 1.0, KE = 1.0)
    add_Se!(BG, S(name * "ATP"))
    add_Se!(BG, S(name * "ADP"))
    add_Ce!(BG, S(name * "E"))
    add_Ce!(BG, S(name * "C"))
    add_Ce!(BG, S(name * "MP"))
    add_Ce!(BG, S(name * "M"))
    for i ∈ 1:6
        add_Bond!(BG, S(name * "B$(i)"))
    end
    add_1J!(BG, Dict(
            S(name * "B1") => true,
            S(name * "B2") => false,
            S(name * "ATP") => true,
            S(name * "M") => true
        ), S(name * "J_M"))
    add_1J!(BG, Dict(
            S(name * "B5") => true,
            S(name * "MP") => false,
            S(name * "ADP") => false,
            S(name * "B6") => false
        ), S(name * "J_MP"))
    add_0J!(BG, Dict(
            S(name * "B3") => true,
            S(name * "B4") => false,
            S(name * "C") => false
        ), S(name * "J_C"))
    add_0J!(BG, Dict(
            S(name * "B6") => true,
            S(name * "E") => false,
            S(name * "B1") => false
        ), S(name * "J_E"))
    add_Re!(BG, S(name * "B2"), S(name * "B3"), S(name * "Re1"))
    add_Re!(BG, S(name * "B4"), S(name * "B5"), S(name * "Re2"))
    kinase_parameters!(BG, name, KM = KM, KE = KE)
end
## Set Kinase Parameters
function kinase_parameters!(BG, name; KM = 1.0, KE = 1.0)
    ATP_potential = affinity_ATP_hydrolysis
    ADP_potential = 0.0
    a = 1000.0
    d = 150.0
    k = 150.0
    KC = exp(ATP_potential) * d / a * KM * KE
    r1 = d / KC
    r2 = k / KC
    BG[S(name * "E")].k = KE
    BG[S(name * "C")].k = KC
    BG[S(name * "M")].k = KM
    BG[S(name * "MP")].k = 0.0
    BG[S(name * "ATP")].Se = ATP_potential
    BG[S(name * "ADP")].Se = ADP_potential
    BG[S(name * "Re1")].r = r1
    BG[S(name * "Re2")].r = r2
end