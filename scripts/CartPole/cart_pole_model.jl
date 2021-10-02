## Setup Project
using DrWatson
@quickactivate "BondGraphModeling"
using Pkg
Pkg.activate(".")
Pkg.instantiate()
## Imports
using BondGraphs
using Symbolics
using ModelingToolkit
using DifferentialEquations
using JLD2
println("Done imports")
## Setup
@variables t e(t) f(t) q(t) p(t)
## Create BondGraph
cart_pole = BondGraph(t)
## Create Elements
add_Se!(cart_pole, :in) 
add_I!(cart_pole, :mc)
add_Bond!(cart_pole, :b3)
add_Bond!(cart_pole, :b4)
add_I!(cart_pole, :mpx; causality = true)
add_Bond!(cart_pole, :b6)
add_Bond!(cart_pole, :b7)
add_I!(cart_pole, :J)
add_Bond!(cart_pole, :b9)
add_Bond!(cart_pole, :b10)
add_Bond!(cart_pole, :b11)
add_Se!(cart_pole, :mpg)
add_I!(cart_pole, :mpy; causality = true)
println("Done Creating Elements")
## Add Junctions
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
    :mpx => false
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
    :mpg => false,
    :mpy => false
    ]), :v_p_y)
println("Done Creating Junctions")
## Add Transformers
@parameters l
@variables θ(t) x(t)
add_MTF!(cart_pole, cos(θ) * l, :b7, :b6, :v_x)
add_MTF!(cart_pole, -sin(θ) * l, :b9, :b10, :v_y)
println("Done Creating Transformers")
## Create Model
generate_model!(cart_pole)
println("Done Creating Model")
## Add Equations for getting the Displacements 
D = Differential(t)
eq1 = D(x) ~ cart_pole[:mc].f
eq2 = D(θ) ~ cart_pole[:J].f
eqns = [eq1, eq2]
sts = [x; θ]
ps = [l]
cart_pole.model = extend(cart_pole.model, ODESystem(eqns, t, sts, ps))
println("Done Updating Equations")
##
resolve_derivative_causality!(cart_pole)
simplify_model!(cart_pole)
##
u0 = [
    x => 0.0,
    θ => randn() * 3.14159 / 180,
    cart_pole[:mc].p => 0.0,
    cart_pole[:J].p => 0.0
    ]
p = [
    l => 1.0,
    cart_pole[:mc].I => 1.0,
    cart_pole[:mpx].I => 0.1,
    cart_pole[:mpy].I => 0.1,
    cart_pole[:mpg].Se => 0.1 * 9.81,
    cart_pole[:in].Se => 0.0, 
    cart_pole[:J].I  => 1.0^2 * 0.1,
    ]
tspan = (0.0, 10.0)
prob = ODAEProblem(cart_pole.model, u0, tspan, p) 
sol = solve(prob, Rodas4())   