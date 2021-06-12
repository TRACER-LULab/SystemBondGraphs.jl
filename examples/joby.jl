## Import Packages
using Symbolics
using DifferentialEquations
using Plots
using ModelingToolkit
## Create Basic Components
@parameters t
## Create a General Bond
function bond(;name)
    @variables e(t) f(t)
    ODESystem(Equation[], t, [e, f], [], defaults=[e => 0.0, f => 0.0], name=name)
end
## All One-Port Elements
function R_element(;name, R::Float64=1.0)
    val = R
    @variables e(t) f(t)
    @parameters R
    eqns = [0.0 ~ R * f - e]
    ODESystem(eqns, t, [e, f], [R], defaults=Dict(R => val, e => 0.0, f => 0.0), name=name)
end

function C_element(;name, C::Float64=1.0)
    val = C
    @variables e(t) f(t) q(t)
    @parameters C
    D = Differential(t)
    eqns = [
            D(q) ~ f,
            0.0 ~ q-C*e
            ]
    ODESystem(eqns, t, [e, f, q], [C], defaults=Dict(C => val, e => 0.0, f => 0.0, q => 0.0), name=name)
end

function I_element(;name, I::Float64=1.0)
    val = I
    @variables e(t) f(t) p(t)
    @parameters I
    D = Differential(t)
    eqns = [
            D(p) ~ e
            0.0 ~ p-f*I
            ]
    ODESystem(eqns, t, [e, f, p], [I], defaults=Dict(I => val, e => 0.0, f => 0.0, p => 0.0), name=name)
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
    [
        0 ~ elems[1].e * (-1)^out[1] - m * elems[2].e * (-1)^out[2], 
        0 ~ m * elems[1].f * (-1)^out[1] - elems[2].f * (-1)^out[2]
    ]
end
function gyrator(;r::Float64, elems::Vector{ODESystem}, out::Vector{Bool})
    [
        0 ~ elems[1].e * (-1)^out[1] - r * elems[2].f * (-1)^out[2],
        0 ~ r * elems[1].f * (-1)^out[1] - elems[2].e * (-1)^out[2]
    ]
end

###################################################################################
###################################################################################
## Elements
# TmassI = Qin*k/b_mass + TmassI
@named Sf1  = flowSource(Sf = Qin)
@named C2   = C_element(C = b_mass)
@named b3   = bond()
@named b4   = bond()
@named C5   = C_element(C = rk_peltier)
@named R6   = R_element(R = rk_plate)
@named b7   = bond()
@named Sf8  = flowSource(Sf = 0.0*t)
@named R9   = R_element(R = rk_plate)
@named b10  = bond()
@named R11  = R_element(R = rk_steak)
# @named Sf12 = flowSource(Sf = sin(3t))
@named b13  = bond()
@named Se14 = effortSource(Se = 72.0)
elements = [Sf1, C2, b3, b4, C5, R6, b7, Sf8, R9, b10, R11, b13, Se14]
## Junctions
J11 = oneJunction(elems = [Sf1, C2, b3], out = [false, true, true])
J01 = zeroJunction(elems = [b3, b4, b7], out = [false, true, true])
J12 = oneJunction(elems = [b4, C5, R6], out = [false, true, true])
J13 = oneJunction(elems = [b7, Sf8, R9, b10], out = [false, false, true, true])
J02 = zeroJunction(elems = [b10, R11, b13], out = [false, true, false, true])
TF = transformer(m = m_conv, elems = [b13, Se14], out = [false, true])
## Equations
 eqns = cat(J11, J01, J12, J13, J02, TF, dims = 1)    
@named sys = ODESystem(eqns, t, systems=elements)
sys = structural_simplify(sys)
## Zero Initial conditions
# u0 = states(sys) |> length |> x->x/2 |> floor |> Int |> zeros
u0 = states(sys) |> length |> zeros
# Solve
prob = ODEProblem(sys, u0, (0, 10.0), sparse = true, jac = true)
sol = solve(prob)
plot(sol, denseplot = false)

## Import Joby Number
k_plate     = 137;
rho_mass    = 135.469; 
v_mass      = 0.0926;
m_mass      = rho_mass*v_mass;
cp_mass     = 0.21;
epsilon     = 0.09; 
k_pyrex     = 0.581;
k_insul     = 0.0231;
k_bot       = k_pyrex + k_insul;
moz         = 8;
cp_steak    = 0.66;
k_steak     = 0.248;
TsurfI      = 60;     
TsprevI     = TsurfI;
TmassI      = TsurfI;
Tnot        = 60;
TinfI       = 72;     
TsteakI     = 32;     
TtprevI     = TsteakI;
Qsun        = 0.1209*3600;
Afl         = 1.2917;     
T           = 0.92;       
b_mass      = m_mass*cp_mass; 
mlb         = moz/16;                 
b_steak     = mlb*cp_steak;           
A_steak     = 0.5;                    
x_steak     = 1/12;                   
rk_steak    = k_steak*A_steak/x_steak;
A_plate     = 0.5625;                 
x_plate     = 0.25/12;                
rk_plate    = k_plate*A_plate/x_plate;
UwindI      = 5.2;
MphToMs     = 0.44704;        
UwindK      = UwindI.*MphToMs;
WToBTUs     = 1/1055.055852;  
hc          = 1.16.*(10.45 - UwindK + 10.0*sqrt(UwindK)).*0.000048919;
m_conv      = hc*A_plate;
sb          = 1.714e-9;          
m_rad       = epsilon*sb*A_plate;
dx          = 0.25/12;         
T52         = 302;             
rk_peltier  = k_bot*A_plate/dx;
Qin         = T*Qsun*Afl; 
time        = 2;      
n           = 100;    
k           = time/n; 
Tpeltier = 60;