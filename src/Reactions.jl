## Add Chemostats - BioModeling
export add_Chemostat!,
        add_Flowstat!,
        add_Ce!,
        add_Re!
"""
Create a linear chemostat
"""
add_Chemostat!(bg, name) = add_Se!(bg, name)
"""
Create a nonliner chemostat
"""
add_Chemostat!(bg, Se, params::Vector{}, name) = add_Se!(bg, Se, params::Vector{}, name)

## Add a Flowstat
"""
Create a linear Flowstat
"""
add_Flowstat!(bg, name) = add_Sf!(bg, name)
"""
Create a nonlinear Flowstat
"""
add_Flowstat!(bg, Se, params::Vector{}, name) = add_Sf!(bg, Se, params::Vector{}, name)

## Concentration Element?
"""
Create a "Linear" Concentration-Element for analysis. Setting Causality to true represents the elements being in derivative causality.
"""
function add_Ce!(bg, name)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv) q(bg.graph_data.iv)
    @parameters k
    D = Differential(bg.graph_data.iv)
    eqns = [
            D(q) ~ f,
            e ~ bg.graph_data.R*bg.graph_data.T*log(k*q)
            ]
    model = ODESystem(
            eqns,
            bg.graph_data.iv,
            [e, f, q],
            [k],
            name = name)
    type = :Ce
    bg[name] = BondGraphNode(model, type)
    nothing
end

"""
Add a reaction component between species "in" and "out" with parameters R, T, r

"""
function add_Re!(bg, name)
    @variables e_in(bg.graph_data.iv) f_in(bg.graph_data.iv)
    @variables e_out(bg.graph_data.iv) f_out(bg.graph_data.iv)
    @parameters r

    eqns = [
        0 ~ f_in + f_out,
        0 ~ f_out - r * (exp(e_in / bg.graph_data.R / bg.graph_data.T) - exp(e_out / bg.graph_data.R / bg.graph_data.T))
    ]

    model = ODESystem(eqns, bg.graph_data.iv, [e_in, e_out, f_in, f_out], [r], name=name)
    type = :Re
    bg[name] = BondGraphNode(model, type)
    nothing
end
