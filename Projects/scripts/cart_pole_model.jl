using DrWatson
@quickactivate "Bond Graph Modeling"
DrWatson.greet()
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
# add_R!(cart_pole, :r1)
add_Bond!(cart_pole, :b3)
add_Bond!(cart_pole, :b4)
add_I!(cart_pole, :mpx; causality=true)
add_Bond!(cart_pole, :b6)
add_Bond!(cart_pole, :b7)
add_I!(cart_pole, :J)
add_Bond!(cart_pole, :b9)
add_Bond!(cart_pole, :b10)
add_Bond!(cart_pole, :b11)
add_Se!(cart_pole, :mpg)
add_I!(cart_pole, :mpy; causality=true)
println("Done Creating Elements")
## Add Junctions
add_1J!(cart_pole, Dict([
    :in => false,
    :mc => true,
    :b3 => false
    # :r1 => true 
    ]), :v_c_x)
add_0J!(cart_pole, Dict([
    :b3 => true,
    :b4 => false,
    :b6 => false
    ]), :J0_x)
add_1J!(cart_pole, Dict([
    :b4 => true,
    :mpx => true
    ]), :v_p_x)
add_1J!(cart_pole, Dict([
    :b7 => true,
    :J => true, 
    :b9 => true
    ]), :v_p_ω)
add_0J!(cart_pole, Dict([
    :b10 => false,
    :b11 => false
    ]), :J0_y)
add_1J!(cart_pole, Dict([
    :b11 => true,
    :mpg => true,
    :mpy => true
    ]), :v_p_y)
println("Done Creating Junctions")
## Add Transformers
@parameters l
@variables θ(t) x(t)
add_MTF!(cart_pole, cos(θ) * l, [θ], [l], [
    :b7 => false,
    :b6 => true
    ], :v_x)
add_MTF!(cart_pole, sin(θ) * l, [θ], [l], [
    :b9 => false,
    :b10 => true
    ], :v_y)
println("Done Creating Transformers")
## Create Model
generate_model!(cart_pole)
println("Done Creating Model")
## Add Equations for getting the Displacements 
eqns = equations(cart_pole.model)
D = Differential(t)
eq1 = D(x) ~ cart_pole.elements[:mc].sys.f
eq2 = D(θ) ~ cart_pole.elements[:J].sys.f
eqns = [eqns; eq1; eq2]
sts = [states(cart_pole.model); x; θ]
ps = [parameters(cart_pole.model); l]
cart_pole.model = ODESystem(eqns, t, sts, ps)
println("Done Updating Equations")
##
resolve_derivative_causality!(cart_pole)
simplify_model!(cart_pole)
println("Done Simplifying Model")
## 
@tagsave datadir("ODEModel", "Cart Pole", "cart_pole_model.jld2")
# save_object("cart_pole.jld2", cart_pole.model)
.