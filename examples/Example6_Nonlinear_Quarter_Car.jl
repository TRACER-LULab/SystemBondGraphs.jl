using SystemBondGraphs
using ModelingToolkit
using OrdinaryDiffEq
using IfElse
using CairoMakie
# Create Bond graph
@parameters t
bg = BondGraph(t)
# Create nonlinear Element Functions
# Parameters
@parameters h, d, U, ks1, ks2, qs0, kt, B
# Tire Spring
Fₜ(e, f, q, t, p) = IfElse.ifelse(q >= 0.0, kt * q, 0)
# Connecting Spring
Fₛ(e, f, q, t, p) = IfElse.ifelse(q <= qs0, ks1 * q, ks1 * qs0 + ks2 * (q - qs0))
# Cubic Damper
Fd(e, f, t, p) = B * f^3
# Velocity Input
vᵢ(e, f, t, p) = IfElse.ifelse(
    U / d * t > 0.0,
    IfElse.ifelse((U / d * t) <= 1, h / d * π * U * cos(π * U / d * t), 0.0),
    0.0,
)
# Create Nonlinear Elements
add_Sf!(bg, vᵢ, [U, d, h], :v_i)
add_C!(bg, :e => Fₜ, [kt], :C_2)
add_C!(bg, :e => Fₛ, [ks1, qs0, ks2], :C_9)
add_R!(bg, :e => Fd, [B], :R_8)
# Create Linear Elements
add_Se!(bg, :mg_us)
add_I!(bg, :m_us)
add_Se!(bg, :mg_s)
add_I!(bg, :m_s)
# Create Junctions
add_0J!(bg, :J01)
add_0J!(bg, :J02)
add_1J!(bg, :J11)
add_1J!(bg, :J12)
add_1J!(bg, :J13)
# Add Connections
add_bond!(bg, :v_i, :J01, :e1)
add_bond!(bg, :J01, :C_2, :e2)
add_bond!(bg, :J01, :J11, :e3)
add_bond!(bg, :J11, :mg_us, :e4)
add_bond!(bg, :J11, :m_us, :e5)
add_bond!(bg, :J11, :J02, :e6)
add_bond!(bg, :J02, :J12, :e7)
add_bond!(bg, :J02, :J13, :e8)
add_bond!(bg, :J12, :C_9, :e9)
add_bond!(bg, :J12, :R_8, :e10)
add_bond!(bg, :J13, :mg_s, :e11)
add_bond!(bg, :J13, :m_s, :e12)

# ## Generate Model
sys = generate_model(bg)
sys, _ = structural_simplify(sys, (inputs(sys), []))

# Add Model Parameters and Initial Conditions Table 5.1
p = Dict{Num,Float64}()
u0 = Dict{Num,Float64}()
fₛ = 1.0
ωₛ = 2 * π * fₛ
p[bg[:m_s].model.I] = 320
p[bg[:m_us].model.I] = p[bg[:m_s].model.I] / 6
p[bg[:v_i].model.U] = 0.9
p[bg[:v_i].model.d] = 1.0
p[bg[:v_i].model.h] = 0.25
p[bg[:C_9].model.ks1] = p[bg[:m_s].model.I] * ωₛ^2
p[bg[:C_9].model.ks2] = 10 * p[bg[:C_9].model.ks1]
p[bg[:C_2].model.kt] = 10 * p[bg[:C_9].model.ks1]
u0[bg[:C_9].model.q] = p[bg[:m_s].model.I] * 9.81 / p[bg[:C_9].model.ks1]
u0[bg[:C_2].model.q] =
    (p[bg[:m_s].model.I] + p[bg[:m_us].model.I]) * 9.81 / p[bg[:C_2].model.kt]
p[bg[:C_9].model.qs0] = 1.3 * u0[bg[:C_9].model.q]
p[bg[:R_8].model.B] = 1500.0
p[bg[:mg_us].model.Se] = 9.81 * p[bg[:m_us].model.I]
p[bg[:mg_s].model.Se] = 9.81 * p[bg[:m_s].model.I]
u0[bg[:m_us].model.p] = 0.0
u0[bg[:m_s].model.p] = 0.0

# ## Create ODAE Problem
tspan = (0.0, 5.0)
prob = ODEProblem(sys, u0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol)
