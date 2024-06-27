
"""
Create a nonlinear R-Element with \$\\Phi_r\$ registered with modelingtoolkit.jl to prevent simplification through the nonlinear function. \$e = \\Phi_r(e, f, t, ps)\$. Params are for any parameters to the nonlinear function.
"""
function add_R!(bg, Φr, ps, name; causality = false)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv)
    if Φr.first == :e
        eqns = [e ~ Φr.second(e, f, bg.graph_data.iv, ps)]
    elseif Φr.first == :f
        eqns = [f ~ Φr.second(e, f, bg.graph_data.iv, ps)]
    end
    model = ODESystem(eqns, bg.graph_data.iv, name = name)
    type = :R
    bg[name] = BondGraphNode(model, type)
    nothing
end

"""
Create a Nonlinear C-Element with \$e = \\phi_c(e, f, q, t, ps)\$. Setting Causality to true represents the elements being in derivative causality.
"""
function add_C!(bg, Φc, ps, name; causality = false)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv) q(bg.graph_data.iv)
    D = Differential(bg.graph_data.iv)
    if Φc.first == :e
        eqns = [
            D(q) ~ f,
            e ~ Φc.second(e, f, q, bg.graph_data.iv, ps), # Integral Causality Form
        ]
    elseif Φc.first == :q
        eqns = [
            D(q) ~ f,
            q ~ Φc.second(e, f, q, bg.graph_data.iv, ps), # Integral Causality Form
        ]
    end
    model = ODESystem(eqns, bg.graph_data.iv, name = name)
    type = :C
    bg[name] = BondGraphNode(model, type)
    nothing
end

"""
Create a Nonlinear I-Element with \$f = \\phi_I(e, f, p, t, ps)\$. Setting Causality to true represents the elements being in derivative causality.
"""
function add_I!(bg, Φi, ps, name; causality = false)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv) p(bg.graph_data.iv)
    D = Differential(bg.graph_data.iv)
    if Φi.first == :f
        eqns = [
            D(p) ~ e,
            f ~ Φi.second(e, f, p, bg.graph_data.iv, ps), # Integral Causality Form
        ]
    elseif Φi.first == :p
        eqns = [
            D(p) ~ e,
            p ~ Φi.second(e, f, p, bg.graph_data.iv, ps), # Integral Causality Form
        ]
    end
    model = ODESystem(eqns, bg.graph_data.iv, [e, f, p], [], name = name)
    type = :I
    bg[name] = BondGraphNode(model, type)
    nothing
end

"""
Create a Nonlinear M-Element with \$p = \\phi_M(e, f, p, q, t, ps)\$. Setting Causality to true represents the elements being in derivative causality.
"""
function add_M!(bg, Φm, ps, name; causality = false)
    @variables e(bg.graph_data.iv) f(bg.graph_data.iv) p(bg.graph_data.iv) q(
        bg.graph_data.iv,
    )
    D = Differential(bg.graph_data.iv)
    if Φm.first == :p
        eqns = [
            D(p) ~ e,
            D(q) ~ f,
            p ~ Φm.second(e, f, p, q, bg.graph_data.iv, ps), # Integral Causality Form
        ]
    elseif Φm.first == :q
        eqns = [
            D(p) ~ e,
            D(q) ~ f,
            q ~ Φm.second(e, f, p, q, bg.graph_data.iv, ps), # Integral Causality Form
        ]
    end
    model = ODESystem(eqns, bg.graph_data.iv, [e, f, p, q], [], name = name)
    type = :M
    bg[name] = BondGraphNode(model, type)
    nothing
end
