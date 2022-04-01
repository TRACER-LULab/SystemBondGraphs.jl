##
using BondGraphs
using ModelingToolkit
using DifferentialEquations
using IfElse
## Create independent variable
@parameters t

## Create Empty Bond Graph
quarter_car = BondGraph(t)

# Create Nonlinear Elements
@parameters h, d, U, ks1, ks2, qs0, kt, B
# Tire Spring
Fₜ(e, f, q, t, p) = IfElse.ifelse(q >= 0.0, kt * q, 0)
# Connecting Spring
Fₛ(e, f, q, t, p) = IfElse.ifelse(q <= qs0, ks1 * q, ks1 * qs0 + ks2 * (q - qs0))
# Cubic Damper
Fd(e, f, t, p) = B * f^3
# Velocity Input
vᵢ(e, f, t, p) = IfElse.ifelse(U / d * t > 0.0,
    IfElse.ifelse(
        (U / d * t) <= 1,
        h / d * π * U * cos(π * U / d * t),
        0.0),
    0.0
)

# Create Nonlinear Elements
add_Sf!(quarter_car, vᵢ, [U, d, h], :v_i)
add_C!(quarter_car, :e => Fₜ, [kt], :C_2)
add_C!(quarter_car, :e => Fₛ, [ks1, qs0, ks2], :C_9)
# add_R!(quarter_car, :e=>Fd, [B], :R_8)

# Create Linear Element
add_R!(quarter_car, :R_8)
add_Se!(quarter_car, :mg_us)
add_I!(quarter_car, :m_us)
add_Se!(quarter_car, :mg_s)
add_I!(quarter_car, :m_s)
add_Se!(quarter_car, :Fin)
# Create Arbitrary Connecting Bonds
add_Bond!(quarter_car, :b3)
add_Bond!(quarter_car, :b6)
add_Bond!(quarter_car, :b7)
add_Bond!(quarter_car, :b10)

# Create Junctions
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
        :R_8 => false,
        :Fin => false),
    :J12)
add_1J!(quarter_car, Dict(
        :b10 => true,
        :mg_s => false,
        :m_s => false),
    :J13)
# Add Information InformationProcessor
@parameters w, b
function IP(input)
    Dict(
        :R_8 =>
            [
                w * (sum(tanh.(input[:C_9])) + sum(tanh.(input[:C_2])) + sum(tanh.(input[:v_i])) + randn()) + b
            ]
    )
end

BondGraphs.add_IP!(quarter_car,
    Dict(
        :C_9 => [quarter_car[:C_9].f],
        :C_2 => [quarter_car[:C_2].f],
        :v_i => [quarter_car[:v_i].f]
    ),
    Dict(
        :R_8 => [quarter_car[:R_8].R]
    ),
    IP,
    :IP
)
## Generate Model
sys = generate_model(quarter_car)

## Simplify Model
sys = structural_simplify(sys)

## Add Model Parameters and Initial Conditions Table 5.1
p = Dict{Num,Any}()
u0 = Dict{Num,Any}()
fₛ = 1.0
ωₛ = 2 * π * fₛ
p[quarter_car[:m_s].I] = 320
p[quarter_car[:m_us].I] = p[quarter_car[:m_s].I] / 6
p[quarter_car[:v_i].U] = 0.9
p[quarter_car[:v_i].d] = 1.0
p[quarter_car[:v_i].h] = 0.25
p[quarter_car[:C_9].ks1] = p[quarter_car[:m_s].I] * ωₛ^2
p[quarter_car[:C_9].ks2] = 10 * p[quarter_car[:C_9].ks1]
p[quarter_car[:C_2].kt] = 10 * p[quarter_car[:C_9].ks1]

u0[quarter_car[:C_9].q] = p[quarter_car[:m_s].I] * 9.81 / p[quarter_car[:C_9].ks1]
u0[quarter_car[:C_2].q] = (p[quarter_car[:m_s].I] + p[quarter_car[:m_us].I]) * 9.81 / p[quarter_car[:C_2].kt]

p[quarter_car[:C_9].qs0] = 1.3 * u0[quarter_car[:C_9].q]
# p[quarter_car[:R_8].B] = 1500.0
p[quarter_car[:mg_us].Se] = 9.81 * p[quarter_car[:m_us].I]
p[quarter_car[:mg_s].Se] = 9.81 * p[quarter_car[:m_s].I]
p[w] = 0.0
p[b] = 1500.0
p[quarter_car[:Fin].Se] = 0.0
u0[quarter_car[:m_us].p] = 0.0
u0[quarter_car[:m_s].p] = 0.0

## Create ODAE Problem 
tspan = (0.0, 4.0)
prob = ODAEProblem(sys, u0, tspan, p)
sol = solve(prob, Rodas4())
##
using Plots
plotly()
plot(sol, vars = [quarter_car[:v_i].f, quarter_car[:C_9].q, quarter_car[:C_2].q, quarter_car[:m_us].f, quarter_car[:m_s].f])
##
using GalacticOptim, Optim, LossFunctions
function loss(x, ps)
    p[w] = x[1]
    p[b] = x[2]
    ps[1] = remake(ps[1], p = ModelingToolkit.varmap_to_vars(p, parameters(sys)))
    sol = solve(ps[1], Rodas4())
    value(L2DistLoss(), ones(length(sol[quarter_car[:C_9].q])) * u0[quarter_car[:C_9].q], sol[quarter_car[:C_9].q], AggMode.Sum()) * 100.0
end
opt_x0 = [1e2, 1e2]
opt_p = [prob]
opt_func = OptimizationFunction(loss, GalacticOptim.AutoFiniteDiff())
opt_prob = OptimizationProblem(opt_func, opt_x0, opt_p)
opt_sol = solve(opt_prob, BFGS())

p[w] = opt_sol.u[1]
p[b] = opt_sol.u[2]
prob = remake(prob, p = ModelingToolkit.varmap_to_vars(p, parameters(sys)))
sol = solve(prob, Rodas4())
plot(sol, vars = [quarter_car[:v_i].f, quarter_car[:C_9].q, quarter_car[:C_2].q, quarter_car[:m_us].f, quarter_car[:m_s].f])
##
R = map(i -> IP(
        Dict(
            :C_9 => sol[quarter_car[:C_9].f][i],
            :C_2 => sol[quarter_car[:C_2].f][i],
            :v_i => sol[quarter_car[:v_i].f][i]
        )
    )[:R_8][1],
    eachindex(sol.t)
)
R = Symbolics.value.(substitute.(R, (Dict(b => p[b], w => p[w]),)))
plot!(sol.t, R, label = "R")