## Phosphorylation Cycle
function phosphorylation_cycle_factory!(BG, name; KM = 1.0, KMP = 1.0, KKin = 1.0)
    # Ports
    add_Ce!(BG, S(name * "M"))
    add_Ce!(BG, S(name * "MP"))
    add_Se!(BG, S(name * "ATP"))
    add_Se!(BG, S(name * "ADP"))
    add_Se!(BG, S(name * "P"))
    add_Ce!(BG, S(name * "Kin"))
    add_Ce!(BG, S(name * "Pho"))
    # Kinase and Phosphatase Reactions
    kinase_factory!(BG, name * "K_", KM = KM, KE = KKin)
    phosphatase_factory!(BG, name * "P_", KMP = KMP)
    phosphorylation_cycle_parameters!(BG, name, KMP = KMP, KM = KM, KKin = KKin)
    # Input and Output ports for cycle
    add_0J!(BG, Dict(
            S(name * "M") => false,
            S(name * "K_M") => false,
            S(name * "P_M") => true
        ), S(name * "J_M"))
    add_0J!(BG, Dict(
            S(name * "MP") => false,
            S(name * "K_MP") => true,
            S(name * "P_MP") => false
        ), S(name * "J_MP"))
    # Swap Species inputs for bond elements
    swap_bond(BG, name * "K_M")
    swap_bond(BG, name * "K_MP")
    swap_bond(BG, name * "P_MP")
    swap_bond(BG, name * "P_M")
    # Swap Components to create Cycle Ports
    swap(BG, name * "K_ATP", name * "ATP")
    swap(BG, name * "K_ADP", name * "ADP")
    swap(BG, name * "K_E", name * "Kin")
    swap(BG, name * "P_P", name * "P")
    swap(BG, name * "P_E", name * "Pho")
end
## Phosphorylation Cycle
function phosphorylation_cycle_parameters!(BG, name; KMP = 1.0, KM = 1.0, KKin = 1.0)
    BG[S(name * "M")].k = KM
    BG[S(name * "MP")].k = KMP
    BG[S(name * "Kin")].k = KKin
    swap_defaults(BG, name, "K_ATP", "ATP")
    swap_defaults(BG, name, "K_ADP", "ADP")
    swap_defaults(BG, name, "K_E", "Kin")
    swap_defaults(BG, name, "P_P", "P")
    swap_defaults(BG, name, "P_E", "Pho")
end