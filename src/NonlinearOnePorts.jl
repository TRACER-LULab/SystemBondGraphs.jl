
"""
Create a nonlinear R-Element with \$\\Phi_r\$ registered with modelingtoolkit.jl to prevent simplification through the nonlinear function. \$e = \\Phi_r(e, f, t, ps)\$. Params are for any parameters to the nonlinear function.
"""
function add_R!(BG::AbstractBondGraph, Φr, ps, name; causality = false)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    if Φr.first == :e
        eqns = [e ~ Φr.second(e, f, BG.model.iv, ps)]
    elseif Φr.first == :f
        eqns = [f ~ Φr.second(e, f, BG.model.iv, ps)]
    end
    sys = ODESystem(eqns, BG.model.iv, name = name)
    props = Dict(
        :type => :R,
        :sys => sys,
        :causality => causality,
        :state_var => []
    )
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

"""
Create a Nonlinear C-Element with \$e = \\phi_c(e, f, q, t, ps)\$. Setting Causality to true represents the elements being in derivative causality.
"""
function add_C!(BG::AbstractBondGraph, Φc, ps, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    D = Differential(BG.model.iv)
    if Φc.first == :e
        eqns = [
            D(q) ~ f,
            e ~ Φc.second(e, f, q, BG.model.iv, ps) # Integral Causality Form
        ]
    elseif Φc.first == :q
        eqns = [
            D(q) ~ f,
            q ~ Φc.second(e, f, q, BG.model.iv, ps) # Integral Causality Form
        ]
    end
    sys = ODESystem(eqns, name = name)
    add_vertex!(BG.graph)
    props = Dict(
        :type => :C,
        :sys => sys,
        :causality => causality,
        :state_var => [sys.q]
    )
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

"""
Create a Nonlinear I-Element with \$f = \\phi_I(e, f, p, t, ps)\$. Setting Causality to true represents the elements being in derivative causality.
"""
function add_I!(BG::AbstractBondGraph, Φi, ps, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv)
    D = Differential(BG.model.iv)
    if Φi.first == :f
        eqns = [
            D(p) ~ e,
            f ~ Φi.second(e, f, p, BG.model.iv, ps) # Integral Causality Form
        ]
    elseif Φi.first == :p
        eqns = [
            D(p) ~ e,
            p ~ Φi.second(e, f, p, BG.model.iv, ps) # Integral Causality Form
        ]
    end
    sys = ODESystem(eqns, BG.model.iv, [e, f, p], [], name = name)
    add_vertex!(BG.graph)
    props = Dict(
        :type => :I,
        :sys => sys,
        :causality => causality,
        :state_var => [sys.p]
    )
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

"""
Create a Nonlinear M-Element with \$p = \\phi_M(e, f, p, q, t, ps)\$. Setting Causality to true represents the elements being in derivative causality.
"""
function add_M!(BG::AbstractBondGraph, Φm, ps, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv) q(BG.model.iv)
    D = Differential(BG.model.iv)
    if Φm.first == :p
        eqns = [
            D(p) ~ e,
            D(q) ~ f,
            p ~ Φm.second(e, f, p, q, BG.model.iv, ps) # Integral Causality Form
        ]
    elseif Φm.first == :q
        eqns = [
            D(p) ~ e,
            D(q) ~ f,
            q ~ Φm.second(e, f, p, q, BG.model.iv, ps) # Integral Causality Form
        ]
    end
    sys = ODESystem(eqns, BG.model.iv, [e, f, p, q], [], name = name)
    props = Dict(
        :type => :M,
        :sys => sys,
        :causality => causality,
        :state_var => [sys.p, sys.q]
    )
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end
