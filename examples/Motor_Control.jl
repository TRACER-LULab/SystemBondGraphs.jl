using BondGraphs
using ModelingToolkit
using DifferentialEquations
##
motor, sys = let 
        @variables t
        motor = BondGraph(t)
        # Add Elements
        add_Se!(motor, :ec)
        add_R!(motor, :Rw)
        add_Bond!(motor, :b3)
        add_Bond!(motor, :b4)
        add_C!(motor, :kτ)
        add_Bond!(motor, :b6)
        add_I!(motor, :J)
        add_Se!(motor, :τd)
        # Add Gyrator
        add_GY!(motor, :b3, :b4, :T)
        # Add Junctions
        add_1J!(motor, Dict(
                :ec => true,
                :Rw => false,
                :b3 => false
                ), :J11)
        add_0J!(motor, Dict(
                :b4 => true,
                :kτ => false,
                :b6 => false
                ), :J01)
        add_1J!(motor, Dict(
                :b6 => true,
                :J => false,
                :τd => true
                ), :J12)
        # Generate Modelalg
        sys = generate_model(motor)
        sys = structural_simplify(sys)
        motor, sys
end;
A, B, C, D, x⃗, u⃗, y⃗ = state_space(motor, sys);
## System Parameters
p = Dict{Num, Float64}()
@parameters m, ωₙ, R
p[motor[:Rw].R] = 1
p[motor[:T].r]  = 0.5
p[m]            = 5/2.2
p[R]            = 4.0*0.0254
p[motor[:J].I]  = p[m]*p[R]^2/2.0
p[ωₙ]           = 2π
p[motor[:kτ].C] = 1/(p[motor[:J].I]*p[ωₙ]^2)
## Update State Space Matrices A, B, C, D
A = Symbolics.value.(substitute.(A, (p,)))
B = Symbolics.value.(substitute.(B, (p,)))
C = Symbolics.value.(substitute.(C, (p,)))
D = Symbolics.value.(substitute.(D, (p,)))
## Convert to Control System
# x⃗̇ = Ax⃗ + Bu⃗
# y⃗ = Cx⃗ + Du⃗
# Output Variable = q
#
C = zeros(size(A, 2))
D = zeros(size(B, 2))
C[x⃗[motor[:kτ].q]] = 1.0
ss(A, B, C, D)