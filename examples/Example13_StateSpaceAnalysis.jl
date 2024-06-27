using OrdinaryDiffEq
using ModelingToolkit: inputs
using SystemBondGraphs
using CairoMakie
using Latexify
@variables t
motor = BondGraph(t)
# Add Elements
add_Se!(motor, :ec)
add_R!(motor, :Rw)
# add_Bond!(motor, :b3)
# add_Bond!(motor, :b4)
add_C!(motor, :kτ)
add_R!(motor, :R)
# add_Bond!(motor, :b6)
add_I!(motor, :J)
add_Se!(motor, :τd)
# Add Gyrator
add_GY!(motor, :T)
# Add Junctions
add_1J!(motor, :j11)
add_bond!(motor, :j11, :Rw, :e1)
add_bond!(motor, :ec, :j11, :e2)
add_bond!(motor, :j11, :T, :e3)
add_0J!(motor, :j01)
add_bond!(motor, :T, :j01, :e4)
add_bond!(motor, :j01, :kτ, :e5)
add_1J!(motor, :j12)
add_bond!(motor, :j01, :j12, :e6)
add_bond!(motor, :j12, :J, :e7)
add_bond!(motor, :τd, :j12, :e8)
add_bond!(motor, :j12, :R, :e9)

# Generate Model
sys = generate_model(motor)
(; A, B, C, D), sys = ModelingToolkit.linearize_symbolic(sys, [motor[:ec].model.Se, motor[:τd].model.Se], [motor[:kτ].model.q, motor[:J].model.p])

x⃗ = [motor[:kτ].model.q, motor[:J].model.p]
u⃗ = [motor[:ec].model.Se, motor[:τd].model.Se]
y⃗ = [motor[:kτ].model.q, motor[:J].model.p]

ss_dict = Dict(
    :A => A,
    :B => B,
    :C => C,
    :D => D,
    :x => x⃗,
    :u => u⃗,
    :y => y⃗
)
p = Dict{Num,Float64}()
@parameters m, ωₙ, R
p[motor[:Rw].model.R] = 1.0
p[motor[:T].model.r] = 0.5
p[m] = 5.0 / 2.2
p[R] = 4.0 * 0.0254
p[motor[:J].model.I] = p[m] * p[R]^2 / 2.0
p[ωₙ] = 1.0
# p[motor_bg[:kτ].C] = 1/(10.0^kτval)
p[motor[:kτ].model.C] = 1 / (p[motor[:J].model.I] * (10^(1.0) * 2 * π)^2)
p[motor[:R].model.R] = 1e-3
A = Symbolics.value.(substitute.(A, (p,))) .|> float
B = Symbolics.value.(substitute.(B, (p,))) .|> float
C = [0 1 / p[motor[:J].model.I]] .|> float
D = [0 0] .|> float
motor_ss = ss(A, B, C, D)
motor_tf = tf(motor_ss)
motor_tf = minreal(motor_tf)
# Analysis
mag, phase, w = bode(motor_tf)
# AR - plot
f, ax, p = lines(w, mag[1,1, :], axis = (xscale=log10, ylabel = "Amplitude Ratio", xlabel = "[rad/s]"), label = "Voltage")
lines!(ax, w, mag[1, 2, :], label="Disturbance Torque")
axislegend()
f
#
f, ax, p = lines(w, phase[1,1, :], axis = (xscale=log10, ylabel = "Phase Angle [°]", xlabel = "[rad/s]"), label = "Voltage")
lines!(ax, w, phase[1, 2, :], label="Disturbance Torque")
axislegend()
f
