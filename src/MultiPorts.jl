## Create C- multiport
"""
**EXPERIMENTAL!!** Create a multiport C-element. Currently only works for completely symbolic models. Needs work to be integrated with the generate model function
"""
function add_C_multiport!(BG::AbstractBondGraph, elements, parameters, name; ϕi = (e, q, params, BG) -> [], ϕk = (e, q, params, BG) -> [])
    # Do the usual setup
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    for i ∈ eachindex(elements)
        add_edge!(BG.graph, BG.graph[name, :name], BG.graph[elements[i].first, :name])
    end

    D = Differential(BG.model.iv)
    # Sort Elements 
    𝐪_1j = filter(x -> x.second == true, elements)
    j = length(𝐪_1j)
    @show j
    𝐞_jp1n = filter(x -> x.second == false, elements)
    n = length(elements)
    @show n
    # Repack elements based on (7.20) & (7.21)
    elements = [𝐪_1j;𝐞_jp1n]
    # Create variable vectors 
    𝐪 = Symbolics.variables(:q, 1:length(elements), T = Symbolics.FnType)
    𝐪 = map(x -> 𝐪[x](BG.model.iv), eachindex(𝐪))
    𝐞 = map(i -> BG[elements[i].first].e, eachindex(elements))
    # Create Derivative Relationships for displacement d/dt(qᵢ) = fᵢ
    deriv_eqns = map(i -> D(𝐪[i]) ~ ParentScope(BG[elements[i].first].f), eachindex(elements))
    𝐞_1j = ϕi(𝐞[j + 1:n], 𝐪[1:j],  parameters, BG)
    e_eqns = map(i -> ParentScope(BG[elements[i].first].e) ~ 𝐞_1j[i], 1:j)

    𝐪_jp1n = ϕk(𝐞[j + 1:n], 𝐪[1:j], parameters, BG)
    q_eqns = map(i -> 𝐪[j + i] ~ 𝐪_jp1n[i], 1:n - j)
    eqns = [deriv_eqns; e_eqns; q_eqns]
    eqns = convert(Vector{Equation}, eqns)
    subsystems = map(i -> BG[elements[i].first], eachindex(elements))
    sys = ODESystem(eqns, BG.model.iv, name = name)
    props = Dict(
            :type => :MPC,
            :sys => sys,
            :causality => false,
            :state_var => []
            )
    set_props!(BG.graph, node_index, props)
    nothing
end

## Create I-multiport
"""
**EXPERIMENTAL!!** Create a multiport I-element. Currently only works for completely symbolic models. Needs work to be integrated with the generate model function
"""
function add_I_multiport!(BG::AbstractBondGraph, elements, parameters, name; ϕi = (p, f, params) -> [], ϕk = (p, f, params) -> [])
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
    𝐩 = Symbolics.variables(:p, 1:length(elements), T = Symbolics.FnType)
    𝐩 = map(x -> 𝐩[x](BG.model.iv), eachindex(𝐩))
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
