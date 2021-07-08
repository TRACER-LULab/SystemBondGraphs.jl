using BondGraphs
using ModelingToolkit
using DifferentialEquations
using Plots
using Symbolics
using SymbolicUtils
using Unitful
## TODO: Implement automated DAE construction with mass matrix to use DAE_index lowering to simplify Problem
# Genereate Mass Matrix from Iterating over the equations. 
## ## ## ## ## ## ## ## ## ## ## ## ## ## 
##  Linear Mass Spring Damper Systems  ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## 

@variables t
msd = BondGraph(t)
# Elements
add_R!(msd, :r1)  
add_C!(msd, :c1)
add_I!(msd, :i1)
add_Se!(msd, 1.0, :se)
# Junctions
add_1J!(msd, Dict([
    :r1 => true, 
    :c1 => true, 
    :i1 => true, 
    :se => false]),
    :J1)
# Creating Model
# Set Problem parameters
msd.elements[:r1].sys.R = 0.00001
msd.elements[:c1].sys.C = 1.0
msd.elements[:i1].sys.I = 200
# Set initial conditions
u0 = [
    msd.elements[:c1].sys.q => 0.0,
    msd.elements[:i1].sys.p => 0.0,
    ]
ps = [
    msd.elements[:c1].sys.C => 0.1,
    msd.elements[:r1].sys.R => 0.00001,
    msd.elements[:i1].sys.I => 3.0
]
# Set TimeSpan
tspan = (0.0, 10.0)
generate_model!(msd, ps)
simplify_model!(msd)
prob = ODEProblem(msd.model, u0, tspan, ps=[1.0 ,100, 3.0])
sol = solve(prob, Rodas5())
plot(sol)

## ## ## ## ## ## ## ## ## ## ## ## ## ##   
##  Nonlinear Quarter Car Figure 5.15  ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Create Bond Graph
@variables t e(t) f(t) q(t) p(t)
quarterCar = BondGraph(t)
# Define Problem values
# Bump
U = 0.9
d = 1.0
h = 0.25
# Masses
ms = 320
mus = ms / 6
# Suspension
fs = 1.0
ωs = 2 * π * fs
ks1 = ms * ωs^2
ks2 = 10 * ks1
kt = 10 * ks1
g = 9.81
qs_init = ms * g / ks1
qus_init = (ms + mus) * g / kt
qs0 = 1.3 * qs_init
# Tire
kt = 10 * ks1
# Damper
B = 1500
# Nonlinear Element Equations
# Bump Velocity Input
v(e, f, t) = ((U / d) * t <= 1) ? (h / d) * π * U * cos(π * U / d * t) : 0.0
@register v(e, f, t)
# Suspension Spring
ϕks(e, q, t) = (q <= qs0) ? ks1 * q : ks1 * qs0 + ks2 * (q - qs0)
@register ϕks(e, q, t)
# Tire Spring
ϕkus(e, q, t) = (q >= 0.0) ? q * kt : 0.0
@register ϕkus(e, q, t)
# Damper Coefficient
ϕb(e, f, t) = B * f^3
@register ϕb(e, f, t)

# Create Elements
add_Sf!(quarterCar, v, [], :Vin)
add_C!(quarterCar, ϕkus, [], :kt)
add_Bond!(quarterCar, :b3)
add_Se!(quarterCar, mus * g, :Fus)
add_I!(quarterCar, :mus)
add_Bond!(quarterCar, :b6)
add_Bond!(quarterCar, :b7)
add_R!(quarterCar, ϕb, [], :b)
add_C!(quarterCar, ϕks, [], :ks)
add_Bond!(quarterCar, :b10)
add_Se!(quarterCar, ms * g, :Fs)
add_I!(quarterCar, :ms)
# Set Linear parameters
quarterCar.elements[:mus].sys.I = mus
quarterCar.elements[:ms].sys.I = ms
# Add Junctions
add_0J!(quarterCar, Dict([
    :Vin => false, 
    :kt => true, 
    :b3 => true]),
    :J01)
add_1J!(quarterCar, Dict([
    :b3 => false, 
    :mus => true, 
    :Fus => true,
    :b6 => true]), 
    :J11)
add_0J!(quarterCar, Dict([
    :b6 => false, 
    :b7 => true, 
    :b10 => true]), 
    :J02)
add_1J!(quarterCar, Dict([
    :b7 => false, 
    :ks => true,
    :b => true]), 
    :J12)
add_1J!(quarterCar, Dict([
    :b10 => false,
    :Fs => true,
    :ms => true]),
    :J13)
# Create Model
generate_model!(quarterCar)
# Set initial conditions
u0 = [
    quarterCar.elements[:ks].sys.q  => qs_init,
    quarterCar.elements[:kt].sys.q  => qus_init,
    quarterCar.elements[:mus].sys.p => 0.0,
    quarterCar.elements[:ms].sys.p  => 0.0
    ]
# Set TimeSpan
tspan = (0.0, 2.0)
simplify_model!(quarterCar)
prob = ODEProblem(quarterCar.model, u0, tspan, sparse=true, jac=true)
sol = solve(prob)
Plots.plot(sol, vars=[quarterCar.elements[:ks].sys.e, quarterCar.elements[:kt].sys.e]) 

## ## ## ## ## ## ## ## ## ## ## ##
##  Algebraic Loop Figure 5.10b  ##
## ## ## ## ## ## ## ## ## ## ## ##

@variables t e(t) f(t) q(t) p(t)
fig5_10 = BondGraph(t)
# Add Elements
add_Se!(fig5_10, 1.0, :F)
add_I!(fig5_10, :m)
add_R!(fig5_10, :b1)
add_R!(fig5_10, :b2)
add_C!(fig5_10, :k)
add_Sf!(fig5_10, t, :Sf)
add_Bond!(fig5_10, :b3)
add_Bond!(fig5_10, :b5)
add_Bond!(fig5_10, :b7)
# Set Linear parameters
fig5_10.elements[:m].sys.I = 1.0
fig5_10.elements[:b1].sys.R = 1.0
fig5_10.elements[:b2].sys.R = 1.0
fig5_10.elements[:k].sys.C = 1.0
# Add Junctions
add_1J!(fig5_10, Dict([
    :F => false,
    :m => true,
    :b3 => true
    ]),
    :J10)
add_0J!(fig5_10, Dict([
    :b3 => false,
    :b1 => true,
    :b5 => true
    ]),
    :J01)
add_1J!(fig5_10, Dict([
    :b5 => false,
    :b2 => true,
    :b7 => true   
    ]),
    :J12)
add_0J!(fig5_10, Dict([
    :b7 => false,
    :k => true,
    :Sf => false
    ]),
    :J02)
# Create system
generate_model!(fig5_10)
u0 = [
    fig5_10.elements[:m].sys.p => 0.0,
    fig5_10.elements[:k].sys.q => 0.0
    ]
# Set TimeSpan
tspan = (0.0, 2.0)
simplify_model!(fig5_10)
prob = ODAEProblem(fig5_10.model, u0, tspan, sparse=true, jac=true)
+sol = solve(prob)
plot(sol)

## ## ## ## ## ## ## ## ## ## ## ## ## ##
##  Derivative Causality Figure 5.12b  ##
## ## ## ## ## ## ## ## ## ## ## ## ## ##
@variables t e(t) f(t) q(t) p(t)
fig5_12 = BondGraph(t)
# Add Elements
add_Se!(fig5_12, sin(t), :S_e)
add_R!(fig5_12, :R_2)
add_I!(fig5_12, :I_3)
add_Bond!(fig5_12, :b4)
add_1J!(fig5_12, Dict([
    :S_e => false,
    :R_2 => true,
:I_3 => true,
    :b4 => true
    ]),
    :J11)
add_Bond!(fig5_12, :b5)
@parameters T
add_GY!(fig5_12, 2.0, Dict([
    :b4 => false,
    :b5 => true
    ]),
    :gy)
add_C!(fig5_12, :C_6)
add_Bond!(fig5_12, :b7)
add_0J!(fig5_12, Dict([
    :b5 => false,
:C_6 => true,
    :b7 => true
    ]),
    :J01)
add_I!(fig5_12, :I_8)
add_R!(fig5_12, :R_9)
add_1J!(fig5_12, Dict([
    :b7 => false,
    :I_8 => true,
    :R_9 => true
    ]),
    :J12)
# Create and Solve
generate_model!(fig5_12)
fig5_12.elements[:I_8].sys.I = 1.0
function get_diff(eqn)
    diff_eqns = []
    alg_eqns = []
    for eqn ∈ eqns
        if !(eqn.lhs isa Float64)
            # Check for differential equations
            if SymbolicUtils.operation(eqn.lhs)  isa Differential
                push!(diff_eqns, eqn)
            # Check for equations with state variables
            else 
                push!(alg_eqns, eqn)
            end
        else 
            push!(alg_eqns, eqn)
end
    end
    return diff_eqns, alg_eqns
end

state_vars = reduce(vcat, map(x -> x.state_var, values(fig5_12.elements)))
eqns = equations(fig5_12.model)
diff_eqns, alg_eqns = get_diff(eqns)
imp_eqn = fig5_12.elements[:C_6].sys.q ~ fig5_12.elements[:C_6].sys.e * fig5_12.elements[:C_6].sys.C
imp_eqn1, eqns1, sub_dict = get_implicit(state_vars, imp_eqn, alg_eqns)
# diff_eqns = map(x->expand_derivatives(substitute(x, Dict([imp_eqn1.lhs => imp_eqn1.rhs]))),diff_eqns)
D = Differential(t)
dt_dc_eqn1 = expand_derivatives(D(imp_eqn1.lhs) ~ D(imp_eqn1.rhs))
diff_eqns = map(x -> substitute(x, Dict([dt_dc_eqn1.lhs => dt_dc_eqn1.rhs])), diff_eqns)
sys = ODESystem([alg_eqns; diff_eqns[2:3]; imp_eqn1], fig5_12.model.iv)
sys = flatten(sys)
sys = alias_elimination(sys)
sys = structural_simplify(sys)
sys = tearing(sys)

function rewrite_eqns(eqns)
    for i ∈ eachindex(eqns)
        eqns[i] = eqns[i].lhs ~ Symbolics.solve_for(eqns[i].rhs - eqns[i].lhs ~ 0.0, eqns[i].lhs)
    end
    return eqns
end

eqns = rewrite_eqns(equations(sys))
sys = ODESystem(eqns, fig5_12.model.iv)
sys = structural_simplify(sys)
sys = tearing(sys)
tspan = (0.0, 2.0)
u0 = [
    fig5_12.elements[:I_8].sys.p => 0.0,
    fig5_12.elements[:I_3].sys.p => 1.0,
    ] |> Dict
ps = [
    fig5_12.elements[:R_9].sys.R => 3.0,
    fig5_12.elements[:C_6].sys.C => 5.0,
    fig5_12.elements[:I_3].sys.I => 2.0,
    fig5_12.elements[:I_8].sys.I => 1.0,
    fig5_12.elements[:R_2].sys.R => 4.0
] |> Dict

prob = ODAEProblem(sys, u0, tspan, ps)
sol = solve(prob, Rodas4())
plot(sol)

