using BondGraphs
using ModelingToolkit
using DifferentialEquations
using IfElse
# # Create Bond graph
@parameters t
quarter_car = BondGraph(t, :qc)
# # Create nonlinear Element Functions
# ## Parameters
@parameters h, d, U, ks1, ks2, qs0, kt, B
# ## Tire Spring
Fₜ(e, f, q, t, p) = IfElse.ifelse(q >= 0.0, kt * q, 0)
# ## Connecting Spring
Fₛ(e, f, q, t, p) = IfElse.ifelse(q <= qs0, ks1 * q, ks1 * qs0 + ks2 * (q - qs0))
# ## Cubic Damper
Fd(e, f, t, p) = B * f^3
# ## Velocity Input
vᵢ(e, f, t, p) = IfElse.ifelse(U / d * t > 0.0,
                                IfElse.ifelse(
                                    (U / d * t) <= 1,
                                    h / d * π * U * cos(π * U / d * t),
                                    0.0),
                                0.0
                                )

# ## Create Nonlinear Elements
add_Sf!(quarter_car, vᵢ, [U, d, h], :v_i)
add_C!(quarter_car, :e => Fₜ, [kt], :C_2)
add_C!(quarter_car, :e => Fₛ, [ks1, qs0, ks2], :C_9)
add_R!(quarter_car, :e => Fd, [B], :R_8)
# ## Create Linear Elements
add_Se!(quarter_car, :mg_us)
add_I!(quarter_car, :m_us)
add_Se!(quarter_car, :mg_s)
add_I!(quarter_car, :m_s)
# ## Create Arbitrary Connecting Bonds
add_Bond!(quarter_car, :b3)
add_Bond!(quarter_car, :b6)
add_Bond!(quarter_car, :b7)
add_Bond!(quarter_car, :b10)
# ## Create Junctions
add_0J!(quarter_car, Dict(
    :v_i => true,
    :C_2 => false,
    :b3 => false
    ),
    :J01)
add_1J!(quarter_car, Dict(
    :b3 => true,
    :mg_us => false,
    :m_us => false,
    :b6 => false
    ),
    :J11)
add_0J!(quarter_car, Dict(
    :b6 => true,
    :b7 => false,
    :b10 => false),
    :J02)
add_1J!(quarter_car, Dict(
    :b7 => true,
    :C_9 => false,
    :R_8 => false),
    :J12)
add_1J!(quarter_car, Dict(
    :b10 => true,
    :mg_s => false,
    :m_s => false),
    :J13)

# ## Generate Model
generate_model!(quarter_car)

# ## Simplify Model
quarter_car.model = structural_simplify(quarter_car.model)

# ## Add Model Parameters and Initial Conditions Table 5.1
p = Dict{Num, Float64}()
u0 = Dict{Num, Float64}()
fₛ = 1.0
ωₛ = 2 * π * fₛ
p[quarter_car[:m_s].I]    = 320
p[quarter_car[:m_us].I]   = p[quarter_car[:m_s].I] / 6
p[quarter_car[:v_i].U]    = 0.9
p[quarter_car[:v_i].d]    = 1.0
p[quarter_car[:v_i].h]    = 0.25
p[quarter_car[:C_9].ks1]  = p[quarter_car[:m_s].I] * ωₛ^2
p[quarter_car[:C_9].ks2]  = 10 * p[quarter_car[:C_9].ks1]
p[quarter_car[:C_2].kt]   = 10 * p[quarter_car[:C_9].ks1]
u0[quarter_car[:C_9].q]   = p[quarter_car[:m_s].I] * 9.81 / p[quarter_car[:C_9].ks1]
u0[quarter_car[:C_2].q]   = (p[quarter_car[:m_s].I] + p[quarter_car[:m_us].I]) * 9.81 / p[quarter_car[:C_2].kt]
p[quarter_car[:C_9].qs0]  = 1.3 * u0[quarter_car[:C_9].q]
p[quarter_car[:R_8].B]    = 1500.0
p[quarter_car[:mg_us].Se] = 9.81 * p[quarter_car[:m_us].I]
p[quarter_car[:mg_s].Se]  = 9.81 * p[quarter_car[:m_s].I]
u0[quarter_car[:m_us].p]  = 0.0
u0[quarter_car[:m_s].p]   = 0.0

# ## Create ODAE Problem
tspan = (0.0, 2.0)
prob = ODAEProblem(quarter_car.model, u0, tspan, p)
sol = solve(prob, Rodas4())
plot(sol)
