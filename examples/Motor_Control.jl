using BondGraphs
using ModelingToolkit
using DifferentialEquations
using LinearAlgebra
using ControlSystems
using Plots
# # Create Bond Graphs and Model
@variables t ω(t)
mat = BondGraph(t, :motor)
add_Se!(mat, :Va)
add_I!(mat, :L)
add_R!(mat, :R)
add_Bond!(mat, :b4)
add_Bond!(mat, :b5)
add_GY!(mat, :b4, :b5, :km)
add_I!(mat, :J)
add_Se!(mat, :τ)
add_R!(mat, :kf)
add_1J!(mat, Dict(
        :Va => true,
        :L => false,
        :R => false,
        :b4 => false
    ), :J11)
add_1J!(mat, Dict(
        :b5 => true,
        :J => false,
        :kf => false,
        :τ => true
    ), :J12)
model = generate_model(mat)
model = structural_simplify(model)
# # Create State-Space Representation for Angular Velocity and Transfer Fucntions
A, B, C, D, x⃗, u⃗, y⃗ = state_space(mat, model)
p = Dict(
    mat[:R].R => 2.0,
    mat[:L].I => 0.5,
    mat[:km].r => 0.1,
    mat[:kf].R => 0.2,
    mat[:J].I => 0.02
)
A = Symbolics.value.(substitute.(A, (p,)))
B = Symbolics.value.(substitute.(B, (p,)))
C = [0 1 / p[mat[:J].I]]
D = [0 0]
motor_ss = minreal(ss(A, B, C, D))
motor_tf = tf(motor_ss)
motor_tf = minreal(motor_tf)
motor_ss, motor_tf
# # Feedforward Controller
y_ff = let
    kff_gain = 1 / motor_tf[1, 1](0)[1]
    kff = ss(tf([kff_gain], [1.0]))
    cl_ff = motor_tf * diagm([kff_gain, 1.0])
    input(x, t) = [1.0, -0.1 * (5.0 <= t <= 10.0)]
    dt = 0.09
    t = 0.0:dt:dt*150
    y_ff, t, x, u = lsim(cl_ff, input, t, dt=dt)
    y_ff
end
# # Feedback Controller
y_fb = let
    kfb_gain = 5.0
    C = append(ss(tf([kfb_gain], [1, 0])), tf([1]))
    OL = motor_ss * C
    fb = ss(zeros(2, 2), zeros(2, 1), zeros(2, 2), [1, 0])
    cl_fb = minreal(feedback(OL, fb))
    input(x, t) = [1.0, -0.1 * (5.0 <= t <= 10.0)]
    dt = 0.09
    t = 0.0:dt:dt*150
    y_fb, t, x, u = lsim(cl_fb, input, t, dt=dt)
    y_fb
end
# # LQR Controller
y_lqr = let
    motor_aug = [1; ss(tf([1], [1, 0]))] * motor_ss[1, 1]
    R = float.(I(1)) * 100.0
    Q11 = 20.0
    Q12 = 21.5
    Q13 = 1.0
    Q22 = 16.0
    Q23 = 4.0
    Q33 = 10.0
    Q = [Q11 Q12 Q13
        Q12 Q22 Q23
        Q13 Q23 Q33]
    K_lqr = lqr(motor_aug, Q, R)
    K_lqr = [44.7214 36.9701 15.1780]
    P = ss(motor_ss.A, motor_ss.B, vcat(motor_ss.C, [0 1], [1 0]), vcat(motor_ss.D, [0 0], [0 0]))
    C = ss(K_lqr * append(tf([1], [1, 0]), tf([1]), tf([1])))
    OL = P * append(C, tf([1]))
    CL = feedback(OL, ss(I(3)), U1=1:3, U2=1:3)
    cl_lqr = CL[1, 1:3:4]
    input(x, t) = [1.0, -0.1 * (5.0 <= t <= 10.0)]
    dt = 0.09
    t = 0.0:dt:dt*150
    y_lqr, t, x, u = lsim(cl_lqr, input, t, dt=dt)
    y_lqr
end
# # Plot the Results
let
    dt = 0.09
    t = 0.0:dt:dt*150
    plot(t, [y_ff', y_fb', y_lqr'], label=["Feedforward" "Feedback" "LQR"])
    plot!(t, ones(length(t)), label="ω-setpoint")
    plot!(t, t -> -0.1 * (5 <= t <= 10), label="Disturbance Torque")
    plot!(legend=:topleft, xlabel="Time", ylabel="Amplitude")
end
