using BondGraphs
using ModelingToolkit
using DifferentialEquations
## Problem Independent Variable
@parameters t

## Create Empty Bondgraph
halfcar = BondGraph(t)

## Create non-linear Input Function
@parameters ω α
v_in(e, f, t, p) = p[1] * sin(p[2] * t)

## Add bondgraph elements
add_Sf!(halfcar, v_in, [α, ω], :sf1)
add_C!(halfcar, :c2)
add_Bond!(halfcar, :b3)
add_I!(halfcar, :i4)
add_Se!(halfcar, :se5)
add_Bond!(halfcar, :b6)
add_Bond!(halfcar, :b24)
add_Bond!(halfcar, :b23)
add_Bond!(halfcar, :b7)
add_Bond!(halfcar, :b10)
add_C!(halfcar, :c8)
add_R!(halfcar, :r9)
add_Se!(halfcar, :se11)
add_I!(halfcar, :i12)
add_Bond!(halfcar, :b13)
add_Bond!(halfcar, :b14)
add_Bond!(halfcar, :b17)
add_C!(halfcar, :c15)
add_R!(halfcar, :r16)
add_Se!(halfcar, :se18)
add_I!(halfcar, :i19)
add_C!(halfcar, :c20)
add_Bond!(halfcar, :b21)
add_Bond!(halfcar, :b27)
add_C!(halfcar, :c22)
add_R!(halfcar, :r25)

## Add Junctions to BondGraph
add_0J!(halfcar, Dict([
    :sf1 => false,
    :c2 => true,
    :b3 => true
    ]),
    :J1)
add_1J!(halfcar, Dict([
    :b3 => false,
    :i4 => true,
    :se5 => false,
    :b6 => true,
    :b23 => true,
    :b24 => true
    ]), 
    :J2)
add_0J!(halfcar, Dict([
    :b6 => false,
    :b7 => true,
    :b10 => true
    ]),
    :J3)
add_1J!(halfcar, Dict([
    :b7 => false,
    :c8 => true,
    :r9 => true
    ]),
    :J4)
add_1J!(halfcar, Dict([
    :b10 => false,
    :se11 => true,
    :i12 => true,
    :b13 => false
    ]),
    :J5)
add_0J!(halfcar, Dict([
    :b13 => true,
    :b14 => true,
    :b17 => false
    ]),
    :J6)
add_1J!(halfcar, Dict([
    :b14 => false,
:c15 => true,
    :r16 => true
    ]),
    :J7)
add_1J!(halfcar, Dict([
    :b17 => true,
    :se18 => false,
    :i19 => true,
    :c20 => true,
    :b21 => false,
    :b27 => false
    ]),
    :J8)
add_0J!(halfcar, Dict([
    :b21 => true,
    :c22 => true,
    :b23 => false
    ]),
    :J9)
add_0J!(halfcar, Dict([
    :b24 => false,
    :r25 => true,
    :b27 => true
    ]),
    :J10)
## Set problem parameters
generate_model!(halfcar)
halfcar.model = structural_simplify(halfcar.model)
## Units and Distributions are Automatically Compatible with the bondgraph package
u0 = [
halfcar[:c2].q  => 0.0,
halfcar[:c8].q  => 0.0,
halfcar[:c15].q => 0.0,
halfcar[:c22].q => 0.0,
halfcar[:i4].p  => 0.0,
halfcar[:i12].p => 0.0,
halfcar[:i19].p => 0.0,
halfcar[:c20].q => 0.0,
]
# Set Parameter Dictionary for parameters of interest
ps = [
halfcar[:sf1].α  => 1.0,
halfcar[:sf1].ω  => 20.0,
halfcar[:c2].C   => 1.0,
halfcar[:i4].I   => 1.0,
halfcar[:se5].Se => 9.81,
halfcar[:c8].C   => 1.0,
halfcar[:r9].R   => 2.0,
halfcar[:i12].I  => 1.0,
halfcar[:se11].Se => 9.8,
halfcar[:c15].C  => 1.0,
halfcar[:r16].R  => 1.0,
halfcar[:i19].I  => 2.0,
halfcar[:se18].Se => 9.8,
halfcar[:c20].C  => 1.0,
halfcar[:c22].C  => 1.0,
halfcar[:r25].R  => 1.0,
]
##
tspan = (0.0, 10.0)
prob = ODAEProblem(halfcar.model, u0, tspan, ps)
sol = solve(prob, Rodas4())