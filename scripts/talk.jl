using DrWatson
@quickactivate "BondGraphModeling"
##
using Symbolics
using ModelingToolkit
using Plots
using DifferentialEquations
##
""" 
Use @register to prevent symbolics from tracing further than necessary
Leverage mathematical structure of bondgraphs and sparsity of the systems
build_function(grad, xs; expression = Val{false}) => avoids runtime world age errors in creating symbolic functions
similar creates an empty of array with size and type of argument
"""
##
@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

D(x)
## 
eqs = [
    D(D(x)) ~ σ * (y - x), 
    D(y) ~ x * (ρ - z),
    D(z) ~ x * y - β * z
    ]
## 
sys = ODESystem(eqs)
sys = ode_order_lowering(sys)
##
u0 = [
    D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0
    ]
p = [
    σ => 28.0,
    ρ => 10.0,
    β => 8 / 3
    ]
tspan = (0.0, 100.0)
## 
prob = ODEProblem(sys, u0, tspan, p, jac=true)
sol = solve(prob, Tsit5())
## 
sol[x]
sol(50.0, idxs=y) # continuous function with respect to a variables
## Check composing verus inheritance for junctions and element structures
# Introduction to acausal Modeling
using OrdinaryDiffEq

@parameters(t)
D = Differential(t)

@connector function Pin(;name)
    sts = @variables v(t) = 1.0 i(t) = 1.0
    ODESystem(Equation[], t, sts, [], name=name)
end

function ModelingToolkit.connect(::Type{Pin}, pins...)
    # KCL ∑ i = 0
    eqns = [
        0 ~ sum(p -> p.i, pins)
        ]
    for i in 1:length(pins) - 1 
        push!(eqns, pins[i].v ~ pins[i + 1].v)
    end
end

function Ground(;name)
    @named g = Pin()
    eqns = [0 ~ g.v]
    compose(ODESystem(eqns, t, [], [], name=name), g)
end

function OnePort(;name)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t) = 1.0 i(t) = 1.0
    eqns = [
        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i
    ]
    compose(ODESystem(eqns, t, sts, []; name=name), p, n)
end

function Resistor(;name, R=1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport # get the voltage and current of the oneport without the namespacing
    ps = @parameters R = R
    eqns = [ v ~ i * R ]
    extend(ODESystem(eqns, t, [], ps; name=name), oneport)
end

function Capacitor(;name, C=1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport # get the voltage and current of the oneport without the namespacing
    ps = @parameters C = C
    eqns = [ D(v) ~ i / C]
    extend(ODESystem(eqns, t, [], ps; name=name), oneport)
end

function ConstantVoltage(;name, V=1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V = V
    eqns = [v ~ V]
    extend(ODESystem(eqns, t, [], ps; name=name), oneport)
end

R = 1.0
C = 1.0
V = 1.0

@named resistor = Resistor(R=R)
@named capacitor = Capacitor(C=C)
@named source = ConstantVoltage(V=V)
@named ground = Ground()

rc_eqns = [
            connect(source.p, resistor.p),
            connect(resistor.n, capacitor.p),
            connect(capacitor.n, source.n, ground)
          ]
@named rc_model = compose(ODESystem(rq_eqns, t), [resistor, capactior, source, ground])

## connect Galatic optim to the MLtoolkit system for the solvers in Perm changing researhc
## optimization of values in a bondgraph system
