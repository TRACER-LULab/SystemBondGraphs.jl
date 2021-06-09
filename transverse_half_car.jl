using BondGraphs
using ModelingToolkit
using Plots
using DifferentialEquations
using Latexify
using Unitful
## Problem Independent Variable and Power Variables
@variables t e(t) f(t) 
@parameters ω α
## Create Empty Bondgraph
halfcar = BondGraph(t)
## Add bondgraph elements
# Create a function for the input
v_in(t, p) = p[1] * sin(p[2] * t)
@register v_in(t, p)
# Setup the Bonds
add_Sf!(halfcar, v_in, [α, ω], :sf1)
add_C!(halfcar, :c2)
add_Bond!(halfcar, :b3)
add_0J!(halfcar, Dict([
    :sf1 => false,
    :c2 => true,
    :b3 => true
    ]),
    :J1)

add_I!(halfcar, :i4)
add_Se!(halfcar, :se5)
add_Bond!(halfcar, :b6)
add_Bond!(halfcar, :b24)
add_Bond!(halfcar, :b23)
add_1J!(halfcar, Dict([
    :b3 => false,
    :i4 => true,
    :se5 => false,
    :b6 => true,
    :b23 => true,
    :b24 => true
    ]), 
    :J2)
add_Bond!(halfcar, :b7)
add_Bond!(halfcar, :b10)
add_0J!(halfcar, Dict([
    :b6 => false,
    :b7 => true,
    :b10 => true
    ]),
    :J3)
add_C!(halfcar, :c8)
add_R!(halfcar, :r9)
add_1J!(halfcar, Dict([
    :b7 => false,
    :c8 => true,
    :r9 => true
    ]),
    :J4)
add_Se!(halfcar, :se11)
add_I!(halfcar, :i12)
add_Bond!(halfcar, :b13)
add_1J!(halfcar, Dict([
    :b10 => false,
    :se11 => true,
    :i12 => true,
    :b13 => false
    ]),
    :J5)
add_Bond!(halfcar, :b14)
add_Bond!(halfcar, :b17)
add_0J!(halfcar, Dict([
    :b13 => true,
    :b14 => true,
    :b17 => false
    ]),
    :J6)
add_C!(halfcar, :c15)
add_R!(halfcar, :r16)
add_1J!(halfcar, Dict([
    :b14 => false,
:c15 => true,
    :r16 => true
    ]),
    :J7)
add_Se!(halfcar, :se18)
add_I!(halfcar, :i19)
add_C!(halfcar, :c20)
add_Bond!(halfcar, :b21)
add_Bond!(halfcar, :b27)
add_1J!(halfcar, Dict([
    :b17 => true,
    :se18 => false,
    :i19 => true,
    :c20 => true,
    :b21 => false,
    :b27 => false
    ]),
    :J8)
add_C!(halfcar, :c22)
add_0J!(halfcar, Dict([
    :b21 => true,
    :c22 => true,
    :b23 => false
    ]),
    :J9)
add_R!(halfcar, :r25)
add_0J!(halfcar, Dict([
    :b24 => false,
    :r25 => true,
    :b27 => true
    ]),
    :J10)
# Set problem parameters
generate_model!(halfcar)
simplify_model!(halfcar)
# Units and Distributions are Automatically Compatible with the bondgraph package
u0 = [
halfcar.elements[:c2].sys.q  => 0.0,
halfcar.elements[:c8].sys.q  => 0.0,
halfcar.elements[:c15].sys.q => 0.0,
halfcar.elements[:c22].sys.q => 0.0,
halfcar.elements[:i4].sys.p  => 0.0,
halfcar.elements[:i12].sys.p => 0.0,
halfcar.elements[:i19].sys.p => 0.0,
halfcar.elements[:c20].sys.q => 0.0,
]
# Set Parameter Dictionary for parameters of interest
ps = [
halfcar.elements[:sf1].sys.α  => 0.01,
halfcar.elements[:c2].sys.C   => 1.0,
halfcar.elements[:i4].sys.I   => 1.0,
halfcar.elements[:se5].sys.Se => 9.81,
halfcar.elements[:c8].sys.C   => 1.0,
halfcar.elements[:r9].sys.R   => 1.0,
halfcar.elements[:i12].sys.I  => 1.0,
halfcar.elements[:se11].sys.Se => 9.8,
halfcar.elements[:c15].sys.C  => 1.0,
halfcar.elements[:r16].sys.R  => 1.0,
halfcar.elements[:i19].sys.I  => 2.0,
halfcar.elements[:se18].sys.Se => 9.8,
halfcar.elements[:c20].sys.C  => 1.0,
halfcar.elements[:c22].sys.C  => 1.0,
halfcar.elements[:r25].sys.R  => 1.0,
halfcar.elements[:sf1].sys.ω  => 20.0
]
tspan = (0.0, 10) 
prob = ODEProblem(halfcar.model, u0, tspan, ps)
sol = solve(prob, Tsit5())
plot(sol)