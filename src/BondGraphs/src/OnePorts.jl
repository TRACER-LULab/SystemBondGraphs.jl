"""

Create an Bond in a BondGraph to connect Junction->Junction, OnePort->TwoPort. It creates an empty bond to be a connector between elements. 

"""
function add_Bond!(BG::BondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    sys = ODESystem(Equation[], BG.model.iv, [e, f], [], name = name)
    # BG.elements[name] = Element(:B, sys, [], false)   
    add_vertex!(BG.graph)
    props = Dict(
                :type => :B,
                :sys => sys,
                :causality => false,
                :state_var => []
            )   
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

"""

Create a Linear R-Element for the bondgraph. Causality will always be false for an R-element since it does not have a "preferred" causality. 

"""
function add_R!(BG::BondGraph, name; causality = false)
    add_vertex!(BG.graph)
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters R
    eqns = [e ~ R * f ]
    sys  = ODESystem(eqns, BG.model.iv, [e, f], [R], name = name)
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
Create a nonlinear R-Element with \$\\Phi_r\$ registered with modelingtoolkit.jl to prevent simplification through the nonlinear function. \$e = \\Phi_r(e, f, t, ps)\$. Params are for any parameters to the nonlinear function. 
"""
function add_R!(BG::BondGraph, Φr, ps; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [e ~ Φr(e, f, BG.model.iv, ps)]
    sys = ODESystem(eqns, BG.model.iv, [e, f], name = name)
    add_vertex!(BG.graph)
    props = Dict(
            :type => :R,
            :sys => sys,
            :causality => causality,
            :state_var => []
            )
    set_prop!(BG.graph, nv(BG.graph), :name, name)
    set_props!(BG.graph, nv(BG.graph), props)
    nothing
end

"""
Create a Linear C-Element for analysis. Setting Causality to true represents the elements being in derivative causality. 
"""
function add_C!(BG::BondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    @parameters C
    D = Differential(BG.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ q / C
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, q], [C], name = name)
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
Create a Nonlinear C-Element with \$e = \\phi_c(e, q, t, ps)\$. Setting Causality to true represents the elements being in derivative causality. 
"""
function add_C!(BG::BondGraph, Φc, ps, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    D = Differential(BG.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ Φc(e, q, BG.model.iv, ps) # Integral Causality Form
            # 0.0 ~ Φc(e, q, BG.model.iv)
            ]
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
Create a Linear I-Element, setting the casuality to true signifies that the element is in derivative causality
"""
function add_I!(BG::BondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv)
    @parameters I
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            f ~ p / I
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, p], [I], name = name)
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
Create a nonlinear I-element, setting the derivative causality to true signifies that the element is in derivative causality. 
"""
function add_I!(BG::BondGraph, Φi, ps, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv)
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            f ~ Φi(p, f, BG.model.iv, ps) # Integral Causality Form
            ]
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
Create a Linear M-element with causality being set to false. Derivative causality for M-elements is still under development. 
"""
function add_M!(BG::BondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv) q(BG.model.iv)
    @parameters M
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            D(q) ~ f,
            p ~ M * q
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, p, q], [M], name = name)
    add_vertex!(BG.graph)
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

"""
Create a nonlinear M-element and with parameters, \$ps\$. Derivative causailty is still under development. 
"""
function add_M!(BG::BondGraph, Φm, ps, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv) q(BG.model.iv)
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            D(q) ~ f,
            p ~ Φm(p, q, BG.model.iv, ps) # Integral Causality Form
            ]
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