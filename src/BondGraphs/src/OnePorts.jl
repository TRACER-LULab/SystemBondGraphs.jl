## Add Generic Bond to Model
function add_Bond!(BG::BondGraph, name)
    @variables e(BG.model.iv) f(BG.model.iv)
    sys = ODESystem(Equation[], BG.model.iv, [e, f], [], name = name)
    BG.elements[name] = Element(:B, sys, [], false)   
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

## Add R_element to Model
function add_R!(BG::BondGraph, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv)
    @parameters R
    eqns = [e ~ R * f ]
    sys  = ODESystem(eqns, BG.model.iv, [e, f], [R], name = name)
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
function add_R!(BG::BondGraph, Φr, params; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv)
    eqns = [e ~ Φr(e, f, BG.model.iv)]
    sys = ODESystem(eqns, BG.model.iv, [e, f], params, name = name)
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

## Add C-element to Model
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
function add_C!(BG::BondGraph, Φc, params, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) q(BG.model.iv)
    D = Differential(BG.model.iv)
    eqns = [
            D(q) ~ f,
            e ~ Φc(e, q, BG.model.iv) # Integral Causality Form
            # 0.0 ~ Φc(e, q, BG.model.iv)
            ]
    sys = ODESystem(eqns, BG.model.iv, [e, f, q], [], name = name)
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

## Add I-element to model
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
function add_I!(BG::BondGraph, Φi, params, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv)
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            f ~ Φi(p, f, BG.model.iv) # Integral Causality Form
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

## Add M-element to model
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
function add_M!(BG::BondGraph, Φm, params, name; causality = false)
    @variables e(BG.model.iv) f(BG.model.iv) p(BG.model.iv) q(BG.model.iv)
    D = Differential(BG.model.iv)
    eqns = [
            D(p) ~ e,
            D(q) ~ f,
            p ~ Φi(p, q, BG.model.iv) # Integral Causality Form
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