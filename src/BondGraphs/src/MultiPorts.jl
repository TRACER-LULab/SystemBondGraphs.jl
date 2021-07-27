## Create C- multiport
function add_C_multiport!(BG::BondGraph, elements, parameters, name; ϕi = (e, q, params) -> [], ϕk = (e, q, params) -> [])
    # Do the usual setup
    D = Differential(BG.model.iv)
    # Sort Elements 
    𝐪_1j = filter(x -> x.second == false, elements)
    j = length(𝐪_1j)
    𝐞_jp1n = filter(x -> x.second == true, elements)
    n = length(elements)
    # Repack elements based on (7.20) & (7.21)
    elements = [𝐪_1j;𝐞_jp1n]
    # Create variable vectors 
    @variables 𝐪[1:length(elements)](BG.model.iv)
    # @variables 𝐟[1:length(elements)](BG.model.iv)
    # @variables 𝐞[1:length(elements)](BG.model.iv)
    𝐞 = map(i -> BG[elements[i].first].e, eachindex(elements))
    # Create Derivative Relationships for displacement d/dt(q_i) = f_i
    deriv_eqns = map(i -> D(𝐪[i]) ~ BG[elements[i].first].f, eachindex(elements))
    # Create Relationships for (7.20) e_i = ϕ_i(q_1j, e_jn, p)
    𝐞_1j = ϕi(𝐞[j + 1:n], 𝐪[1:j],  parameters)
    e_eqns = map(i -> BG[elements[i].first].e ~ 𝐞_1j[i], 1:j)
    𝐪_jp1n = ϕk(𝐞[j + 1:n], 𝐪[1:j], parameters)
    q_eqns = map(i -> 𝐪[j + i] ~ 𝐪_jp1n[i], 1:n - j)
    eqns = [deriv_eqns; e_eqns; q_eqns]
    eqns = convert(Vector{Equation}, eqns)
    subsys = map(i -> BG[elements[i].first], eachindex(elements))
    sys = compose(ODESystem(eqns, BG.model.iv, collect(𝐪), [], name = name), subsys)
    BG.elements[name] = Element(:C, sys, collect(𝐪), false)
    nothing
end

## Create I-multiport
function add_I_multiport!(BG::BondGraph, elements, parameters; ϕi = (p, f, params) -> [], ϕk = (p, f, params) -> [], name)
    # Do the usual setup
    D = Differential(BG.model.iv)
    # Sort Elements 
    𝐩_1j = filter(x -> x.second == false, elements)
    j = length(𝐩_1j)
    𝐟_jp1n = filter(x -> x.second == true, elements)
    n = length(elements)
    # Repack elements based on (7.20) & (7.21)
    elements = [𝐩_1j;𝐟_jp1n]
    # Create variable vectors 
    @variables 𝐩[1:length(elements)](BG.model.iv)
    𝐟 = map(i -> BG[elements[i].first].f, eachindex(elements))
    # Create Derivative Relationships for displacement d/dt(q_i) = f_i
    deriv_eqns = map(i -> D(𝐩[i]) ~ BG[elements[i].first].e, eachindex(elements))
    # Create Relationships for (7.20) e_i = ϕ_i(q_1j, e_jn, p)
    𝐟_1j = ϕi(𝐩[1:j], 𝐟[j + 1:n], parameters)
    e_eqns = map(i -> BG[elements[i].first].f ~ 𝐟_1j[i], 1:j)
    𝐩_jn = ϕk(𝐩[1:j], 𝐟[j + 1:n], parameters)
    q_eqns = map(i -> 𝐩[j + i] ~ 𝐩_jn[i], 1:n - j)
    return [deriv_eqns; e_eqns; q_eqns]
end