## Imports
using BondGraphs
using Symbolics
using ModelingToolkit
using DifferentialEquations
println("Done imports")
## Create BondGraph
@variables t
cart_pole = BondGraph(t, :cart_pole)
# Create Elements
add_Se!(cart_pole, :in)
add_I!(cart_pole, :mc)
add_Bond!(cart_pole, :b3)
add_Bond!(cart_pole, :b4)
add_I!(cart_pole, :mpx)
add_Bond!(cart_pole, :b6)
add_Bond!(cart_pole, :b7)
add_I!(cart_pole, :J)
add_Bond!(cart_pole, :b9)
add_Bond!(cart_pole, :b10)
add_Bond!(cart_pole, :b11)
add_Se!(cart_pole, :mpg)
add_I!(cart_pole, :mpy)
println("Done Creating Elements")
# Add Junctions
add_1J!(cart_pole, Dict([
        :in => true,
        :mc => false,
        :b3 => true
    ]), :v_c_x)
add_0J!(cart_pole, Dict([
        :b3 => false,
        :b4 => true,
        :b6 => true
    ]), :J0_x)
add_1J!(cart_pole, Dict([
        :b4 => false,
        :mpx => true
    ]), :v_p_x)
add_1J!(cart_pole, Dict([
        :b7 => false,
        :J => false,
        :b9 => false
    ]), :v_p_ω)
add_0J!(cart_pole, Dict([
        :b10 => true,
        :b11 => false
    ]), :J0_y)
add_1J!(cart_pole, Dict([
        :b11 => true,
        :mpg => true,
        :mpy => false
    ]), :v_p_y)
@variables θ(t) x(t)
@parameters l
add_MTF!(cart_pole, cos(GlobalScope(θ)) * GlobalScope(l), :b7, :b6, :v_x)
add_MTF!(cart_pole, -sin(GlobalScope(θ)) * GlobalScope(l), :b9, :b10, :v_y)
D = Differential(cart_pole.model.iv)
eqns = [
    D(θ) ~ cart_pole[:J].f,
    D(x) ~ cart_pole[:mc].f,
]
@named theta_model = ODESystem(eqns, cart_pole.model.iv, [θ, x], [])
sys = generate_model(cart_pole)
sys = extend(sys, theta_model)
sys = structural_simplify(sys)
##
# model = BondGraphs.graph_to_model(cart_pole)
# @variables t
# sys = derivative_casuality(cart_pole, skip = [cart_pole[:v_x].θ]);
# @named sys = ODESystem(expand.(ModelingToolkit.get_eqs(sys)), t)
##
# mpx₊p(t) ~ mpx₊I*(mc₊p(t) / mc₊I - l*((-((-l*J₊p(t)*sin(θ(t))) / J₊I)) / (l*sin(θ(t))))*cos(θ(t)))
##
sys = alias_elimination(sys)
sys = initialize_system_structure(sys)
# model = structural_simplify(model)
##
# @variables mc₊p
u0 = [
    x => 0.0,
    θ => randn() * 3.14159 / 180,
    cart_pole[:mpx].f => 0.0,
    cart_pole[:b10].f => 0.0,
]
p = [
    l => 1.0,
    cart_pole[:mc].I => 1.0,
    cart_pole[:mpx].I => 1.0,
    cart_pole[:mpy].I => 1.0,
    cart_pole[:J].I => 1.0^2 * 1.0,
    cart_pole[:mpg].Se => 1.0 * 9.81,
    cart_pole[:in].Se => 0.0,
]
tspan = (0.0, 10.0)
new_prob = ODAEProblem(sys, u0, tspan, p)
sol = solve(new_prob, Tsit5())
plot(sol, vars = [x, θ])

## TESTING
@variables b9₊f(t) mc₊f(t) mc₊e(t) J₊e(t) θ(t) J₊p(t) b6₊e(t) mc₊p(t) mpy₊e(t) in₊e(t) b7₊e(t)
@variables b9₊e(t) b10₊e(t) mpg₊e(t) b6₊f(t) mpx₊f(t) b10₊f(t) mpx₊p(t) mpy₊p(t)
@parameters mpx₊I, l, mc₊I, J₊I, mpy₊I, in₊Se, mpg₊Se
D = Differential(t)
eqs = [
    D(θ) ~ b9₊f
    D(x) ~ mc₊f
    D(mc₊p) ~ mc₊e
    D(J₊p) ~ J₊e
    0 ~ (mpx₊I * ((l * (mc₊e - l * mc₊I * ((l * J₊e * sin(θ) + l * J₊p * b9₊f * cos(θ) - J₊I * (l^2) * ((-((-l * J₊p * sin(θ)) / J₊I)) / ((l^2) * (sin(θ)^2))) * b9₊f * cos(θ) * sin(θ)) / (J₊I * l * sin(θ))) * cos(θ)) * sin(θ) - l * mc₊I * ((-l * J₊p * sin(θ)) / J₊I) * b9₊f * sin(θ)) / (l * sin(θ)))) / mc₊I - b6₊e
    0 ~ (-l * mpy₊I * (mpx₊I * mc₊p - mc₊I * ((l * mpx₊I * J₊p * cos(θ) - J₊I * mpx₊I * (mc₊p / mc₊I)) / (-J₊I))) * b9₊f * cos(θ) - mc₊I * mpx₊I * mpy₊I * (l^2) * ((J₊I * mpx₊I * mc₊e + J₊I * (l^2) * (mc₊I^2) * (mpx₊I^2) * ((mpx₊I * mc₊p - mc₊I * ((l * mpx₊I * J₊p * cos(θ) - J₊I * mpx₊I * (mc₊p / mc₊I)) / (-J₊I))) / ((l^2) * (mc₊I^2) * (mpx₊I^2) * (cos(θ)^2))) * b9₊f * cos(θ) * sin(θ) - mc₊I * ((mc₊I * (l * mpx₊I * J₊p * b9₊f * sin(θ) - l * mpx₊I * J₊e * cos(θ)) + J₊I * mpx₊I * mc₊e) / mc₊I)) / (J₊I * l * mc₊I * mpx₊I * cos(θ))) * cos(θ) * sin(θ)) / (l * mc₊I * mpx₊I * cos(θ)) - mpy₊e
    0 ~ b6₊e + in₊e - mc₊e
    0 ~ -J₊e - b7₊e - b9₊e
    0 ~ mpg₊e + mpy₊e - b10₊e
    0 ~ b6₊f + mpx₊f - mc₊f
    0 ~ b7₊e - l * b6₊e * cos(θ)
    0 ~ l * b9₊f * cos(θ) - b6₊f
    0 ~ l * b10₊e * sin(θ) + b9₊e
    0 ~ -b10₊f - l * b9₊f * sin(θ)
    0 ~ in₊Se - in₊e
    0 ~ mc₊p / mc₊I - mc₊f
    0 ~ mpx₊p / mpx₊I - mpx₊f
    0 ~ J₊p / J₊I - b9₊f
    0 ~ mpg₊Se - mpg₊e
    0 ~ mpy₊p / mpy₊I - b10₊f
]

@named sys = ODESystem(eqs, t)
sys = initialize_system_structure(sys)
prob = ODAEProblem(sys, u0, tspan, p)
sol = solve(prob, Rodas4(), abstol = 1e-8, restol = 1e-8)
