using BondGraphs
using ModelingToolkit
using DifferentialEquations
##
dc = let
    @variables t
    dc = BondGraph(t)
    # Add Elements
    add_Se!(dc, :ec)
    add_R!(dc, :Rw)
    add_I!(dc, :L)
    add_Bond!(dc, :b4)
    add_Bond!(dc, :b5)
    @parameters r
    add_MGY!(dc, r, :b4, :b5, :T)
    add_C!(dc, :kτ, causality = true)
    add_Bond!(dc, :b7)
    add_I!(dc, :J)
    add_R!(dc, :bτ)
    # Add Junctions
    add_1J!(dc, Dict(
            :ec => true,
            :Rw => false,
            :L => false,
            :b4 => false
        ), :J11)
    add_0J!(dc, Dict(
            :b5 => true,
            :kτ => false,
            :b7 => false
        ), :J01)
    add_1J!(dc, Dict(
            :b7 => true,
            :J => false,
            :bτ => false
        ), :J12)
    dc
end
##
sys = derivative_casuality(dc)
sys = structural_simplify(sys)
## 
u0 = [
    dc[:L].p => 0.0,
    dc[:J].p => 0.0,
]
p = [
    dc[:T].r => 0.5,
    dc[:ec].Se => 1.0,
    dc[:Rw].R => 1.0,
    dc[:L].I => 1.0,
    dc[:kτ].C => 1.0,
    dc[:J].I => 1.0,
    dc[:bτ].R => 1.0
]
tspan = (0.0, 1.0)
prob = ODAEProblem(sys, u0, tspan, p)
sol = solve(prob, restol = 1e-8)