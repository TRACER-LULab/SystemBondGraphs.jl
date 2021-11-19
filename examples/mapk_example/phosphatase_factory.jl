## Phosphatase Factory
function phosphatase_factory!(BG, name; KMP = 1.0, KE = 1.0)
    add_Se!(BG, S(name * "P"))
    add_Ce!(BG, S(name * "E"))
    add_Ce!(BG, S(name * "MP"))
    add_Ce!(BG, S(name * "M"))
    add_Ce!(BG, S(name * "C"))
    for i ∈ 1:6
        add_Bond!(BG, S(name * "B$(i)"))
    end
    add_1J!(BG, Dict(
            S(name * "B1") => true,
            S(name * "B2") => false,
            S(name * "MP") => true
        ), S(name * "J_MP"))
    add_1J!(BG, Dict(
            S(name * "B5") => true,
            S(name * "B6") => false,
            S(name * "M") => false,
            S(name * "P") => false
        ), S(name * "J_M"))
    add_0J!(BG, Dict(
            S(name * "B3") => true,
            S(name * "C") => false,
            S(name * "B4") => false
        ), S(name * "J_C"))
    add_0J!(BG, Dict(
            S(name * "B6") => true,
            S(name * "B1") => false,
            S(name * "E") => false
        ), S(name * "J_E"))
    add_Re!(BG, S(name * "B2"), S(name * "B3"), S(name * "Re1"))
    add_Re!(BG, S(name * "B4"), S(name * "B5"), S(name * "Re2"))
    phosphatase_parameters!(BG, name, KMP = KMP, KE = KE)
end
## Phosphatase Parameters
function phosphatase_parameters!(BG, name; KMP = 1.0, KE = 1.0)
    P_potential = 0.0
    a = 1000.0
    d = 150.0
    k = 150.0
    KC = d / a * KMP * KE
    r1 = d / KC
    r2 = k / KC
    BG[S(name * "E")].k = KE
    BG[S(name * "C")].k = KC
    BG[S(name * "M")].k = 0.0
    BG[S(name * "MP")].k = KMP
    BG[S(name * "P")].Se = P_potential
    BG[S(name * "Re1")].r = r1
    BG[S(name * "Re2")].r = r2
end