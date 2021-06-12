## Import Packages
using Symbolics
using DifferentialEquations
using Plots
using ModelingToolkit
## Create a General Bond
function bond(;name)
    @variables e(t) f(t)
    ODESystem(Equation[], t, [e, f], [], defaults=[e => 0.0, f => 0.0], name=name)
end
## All One-Port Elements
function R_element(;name, R)
    val = R
    @variables e(t) f(t)
    @parameters R
    eqns = [0.0 ~ R * f - e]
    # eqns = [0.0 ~ ΦR(f) - e]
    ODESystem(eqns, t, [e, f], [R], defaults=Dict(R => val, e => 0.0, f => 0.0), name=name)
end

function C_element(;name, C,)
    val = C
    @variables e(t) f(t) q(t)
    @parameters C
    D = Differential(t)
    eqns = [
            D(q) ~ f,
            0.0 ~ q-C*e
            # q ~ C*e
            ]
    ODESystem(eqns, t, [e, f, q], [C], observed = [0.0 ~ q-C*e], defaults=Dict(C => val, e => 0.0, f => 0.0, q => 0.0), name=name)
end

function I_element(;name, I)
    val = I
    @variables e(t) f(t) p(t)
    @parameters I
    D = Differential(t)
    eqns = [
            D(p) ~ e,
            0.0 ~ p-f*I
            # p ~ f*I
            ]
    ODESystem(eqns, t, [e, f, p], [I], observed = [0.0 ~ p-f*I], defaults=Dict(I => val, e => 0.0, f => 0.0, p => 0.0), name=name)
end

function effortSource(;name, Se)
    @variables e(t) f(t)
    eqns = [0.0 ~ e - Se]
    ODESystem(eqns, t, [e, f], [], name=name)
end

function flowSource(;name, Sf)
    @variables e(t) f(t)
    eqns = [0.0 ~ f - Sf]
    ODESystem(eqns, t, [e, f], [], name=name)
end
## Create Two Port Elements
function oneJunction(;elems::Vector{ODESystem}, out::Vector{Bool})
    eqns = [
            0.0 ~ sum(map(i -> elems[i].e * (-1).^(out[i]), eachindex(elems))) # Sum of all efforts is 0
            ]
    for i ∈ 1:length(elems) - 1
        push!(eqns, 0.0 ~ elems[i].f - elems[i + 1].f) # flow equality
    end
    return eqns
end
## Zero Junction
function zeroJunction(;elems::Vector{ODESystem}, out::Vector{Bool})
    eqns = [
            0.0 ~ sum(map(x -> elems[x].f * (-1).^(out[x]), eachindex(elems))) # Sum of all flows is 0
            ]
    for i ∈ 1:length(elems) - 1
        push!(eqns, 0.0 ~ elems[i].e - elems[i + 1].e) # effort equality
    end
    return eqns
end
## Transformer and Gyrator
function transformer(;m::Float64, elems::Vector{ODESystem}, out::Vector{Bool})
    eqns = [
        0 ~ elems[1].e * (-1)^out[1] - m * elems[2].e * (-1)^out[2], 
        0 ~ m * elems[1].f * (-1)^out[1] - elems[2].f * (-1)^out[2]
    ]
end
function gyrator(;r::Float64, elems::Vector{ODESystem}, out::Vector{Bool})
    eqns = [
        0 ~ elems[1].e * (-1)^out[1] - r * elems[2].f * (-1)^out[2],
        0 ~ r * elems[1].f * (-1)^out[1] - elems[2].e * (-1)^out[2]
    ]
end

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## selection of bond graphs from the book ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## Mass Spring Damper systems
# Elements
@named R1 = R_element(R=1.0)
@named C1 = C_element(C=1 / 20.0)
@named I1 = I_element(I=0.2)
@named Se = effortSource(Se=sin(t))
# Junctions
J1 = oneJunction(elems=[R1, C1, I1, Se], out=[true, true, true, false])
# Equations
eqns = cat(J1, dims=1)
@named sys = ODESystem(eqns, t, systems=[R1, C1, I1, Se])
# Simplify
sys = structural_simplify(sys)
# Zero Initial Conditions
u0 = states(sys) |> length |> x->x/2 |> floor |> Int |> zeros
# Solve
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
plot(sol, denseplot=true)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## Example From Figure 6.15 from Textbook page 251
# Elements
@named R1 = R_element(R=1.0)
@named C1 = C_element(C=2.0)
@named I1 = I_element(I=3.0)
@named Sf = flowSource(Sf=sin(t))
@named Se = effortSource(Se=cos(2t))
@named bond3 = bond()
# Junctions
J0_eqns = zeroJunction(elems=[Sf, I1, bond3], out=[false, true, true])
J1_eqns = oneJunction(elems=[bond3, C1, R1, Se], out=[false, true, true, true])
# Equations
eqns = cat(J0_eqns, J1_eqns, dims=1)
@named sys = ODESystem(eqns, t, systems=[R1, C1, I1, Sf, Se, bond3])
# Simplify
sys = structural_simplify(sys)
# Zero Initial Conditions
u0 = states(sys) |> length |> x->x/2 |> floor |> Int |> zeros
# Solve
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
plot(sol, denseplot=false)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## Quater Car Model from Figure 6.47 Page 283
# Elements
@named bs = R_element(R=1.0) 
@named kt = C_element(C=1.0)
@named ks = C_element(C=1.0)
@named mu = I_element(I=1.0)
@named ms = I_element(I=1.0)
@named Fc = effortSource(Se=sin(t))
@named vin = flowSource(Sf=cos(5t))
@named b3 = bond()
@named b5 = bond()
@named b7 = bond()
# Junctions
J01 = zeroJunction(elems=[vin, kt, b3], out=[false, true, true])
J11 = oneJunction(elems=[b3, mu, b5], out=[false, true, true])
J02 = zeroJunction(elems=[b5, b7, ms], out=[false, true, true])
J12 = oneJunction(elems=[b7, Fc, ks, bs], out=[false, true, true, true])
# Equations
eqns = cat(J01, J11, J02, J12, dims=1)
@named sys = ODESystem(eqns, t, systems=[bs, kt, ks, mu, ms, Fc, vin, b3, b5, b7])
# Simplify
sys = structural_simplify(sys)
# Zero Initial Conditions
u0 = states(sys) |> length |> x->x/2 |> floor |> Int |> zeros
# Solve
prob = ODAEProblem(sys, u0, (0, 20.0))
sol = solve(prob, Rodas4())
plot(sol)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Figure 6.52 On Page 290
# Elements
Rw = 1.0
T = 0.5
J = 0.0117
@named Se = effortSource(Se = sin((T^2/Rw/J)*t))
@named Rw = R_element(R = Rw)
@named J = I_element(I = J)
@named τd = effortSource(Se = 0t)
@named bond2 = bond()
@named bond3 = bond()
# Junctions
Junc11 = oneJunction(elems = [Se, Rw, bond2], out = [false, true, true])
Gy = gyrator(r = T, elems = [bond2, bond3], out = [false, true])
Junc12 = oneJunction(elems = [bond3, J, τd], out = [false, true, true])
# Equations
eqns = cat(Junc11, Gy, Junc12, dims = 1)
@named sys = ODESystem(eqns, t, systems = [Se, Rw, J, τd, bond2, bond3])
# Simplify
sys = structural_simplify(sys)
# Zero Initial Conditions
u0 = states(sys) |> length |> x->x/2 |> floor |> Int |> zeros
# Solve
prob = ODAEProblem(sys, u0, (0, 1.0))
sol = solve(prob, Rodas4())
plot(sol)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Figure 6.57 from Page 293
# Elements
Rw = 1.0
T = 0.5
kτ = 1/10.0
J = 0.0117
@named ec = effortSource(Se = sin(T^2/Rw/J*t))
@named Rw = R_element(R = Rw)
@named kτ = C_element(C = kτ)
@named J = I_element(I=J)
@named τd = effortSource(Se = sin(10t))
@named bond3 = bond()
@named bond4 = bond()
@named bond6 = bond()
# Junctions 
Junc11 = oneJunction(elems = [ec, Rw, bond3], out = [false, true, true])
GY = gyrator(r = T, elems = [bond3, bond4], out = [false, true])
Junc01 = zeroJunction(elems = [bond4, kτ, bond6], out = [false, true, true])
Junc12 = oneJunction(elems = [bond6, J, τd], out = [false, true, false])
# Equations
eqns = cat(Junc11, GY, Junc01, Junc12, dims = 1)
@named sys = ODESystem(eqns, t, systems = [ec, Rw, kτ, J, τd, bond3, bond4, bond6])
# Simplify
sys = structural_simplify(sys)
# Zero Initial Conditions
u0 = states(sys) |> length |> x->x/2 |> floor |> Int |> zeros
# Solve
prob = ODAEProblem(sys, u0, (0, 1.0))
sol = solve(prob, Rodas4())
plot(sol, vars = [J.f])

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Figure 5.10 from Page 184 ~~ Algebraic Loops
# Elements
@named Se = effortSource(Se = sin(t))
@named I2 = I_element(I = 1.0)
@named b3 = bond()
@named R4 = R_element(R = 1.0)
@named b5 = bond()
@named R6 = R_element(R = 1.0)
@named b7 = bond()
@named C8 = C_element(C = 1.0)
@named Sf = flowSource(Sf = sin(10t))
# Junctions
Junc11 = oneJunction(elems = [Se, I2, b3], out = [false, true, true])
Junc01 = zeroJunction(elems = [b3, R4, b5], out = [false, true, true])
Junc12 = oneJunction(elems = [b5, R6, b7], out = [false, true, true])
Junc02  = zeroJunction(elems = [b7, C8, Sf], out = [false, true, false])
# Equations
eqns = cat(Junc11, Junc01, Junc12, Junc02, dims = 1)
@named sys = ODESystem(eqns, t, systems = [Se, I2, b3, R4, b5, R6, b7, C8, Sf])
# Simplify
sys = structural_simplify(sys)
# Zero Initial Conditions
u0 = states(sys) |> length |> x->x/2 |> floor |> Int |> zeros
# Solve
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())
plot(sol, denseplot = false)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Figure 5.12 c from Page 190 ~~ Derivative Casuality
# Elements
@named Se = effortSource(Se = sin(t))
@named R2 = R_element(R = 1.0)
@named I3 = I_element(I = 2.0)
@named b4 = bond()
@named b5 = bond()
@named C6 = C_element(C = 3.0)
@named b7 = bond()
@named I8 = I_element(I = 4.0)
@named R9 = R_element(R = 5.)
# Junctions
Junc11 = oneJunction(elems = [Se, R2, I3, b4], out = [false, true, true, true])
GY = gyrator(r = 0.5, elems = [b4, b5], out = [false, true])
Junc01 = zeroJunction(elems = [b5, C6, b7], out = [false, true, true])
Junc12 = oneJunction(elems = [b7, I8, R9], out = [false, true, true])
# Equations
eqns = cat(Junc11, GY, Junc12, Junc01, dims = 1)
@named sys = ODESystem(eqns, t, systems = [Se, R2, I3, C6, I8, R9, b4, b5, b7])
sys = flatten(sys) 
sys = tearing(sys)
sys = alias_elimination(sys)
sys = structural_simplify(sys)
sys = initialize_system_structure(sys)
# Simplify ~ Not working for Derivative Casuality
# Zero Initial Conditions
u0 = states(sys) |> length  |> zeros
# Solve
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob)
plot(sol)