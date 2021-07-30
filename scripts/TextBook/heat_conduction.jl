using Symbolics: getindex
using BondGraphs
using ModelingToolkit
using DifferentialEquations
using Plots
using Symbolics
using SymbolicUtils
using Unitful
using UnitfulRecipes
##
@variables t
wall = BondGraph(t)
##
add_Se!(wall, :T_i)
add_R!(wall, :R_si)
add_Bond!(wall, :B_3)
add_C!(wall, :C_l)
add_Bond!(wall, :B_5)
add_R!(wall, :R_l)
add_Bond!(wall, :B_7)
add_C!(wall, :C_ll)
add_Bond!(wall, :B_9)
add_R!(wall, :R_ll)
add_Bond!(wall, :B_11)
add_C!(wall, :C_lll)
add_Bond!(wall, :B_13)
add_R!(wall, :R_se)
add_Se!(wall, :T_e)
## Add Junctions
add_1J!(wall, Dict([
    :T_i => false,
    :R_si => true,
    :B_3 => true
    ]),
    :J1_1)
add_0J!(wall, Dict([
    :B_3 => false,
    :C_l => true,
    :B_5 => true
    ]), 
    :J0_1)
add_1J!(wall, Dict([
    :B_5 => false,
    :R_l => true,
    :B_7 => true
    ]),
    :J1_2)
add_0J!(wall, Dict([
    :B_7 => false,
    :C_ll => true,
    :B_9 => true
    ]),
    :J0_2)
add_1J!(wall, Dict([
    :B_9 => false,
    :R_ll => true,
    :B_11 => true  
    ]),
    :J1_3)
add_0J!(wall, Dict([
     :B_11 => false,
     :C_lll => true,
     :B_13 => false
    ]),
    :J0_3)
add_1J!(wall, Dict([
    :B_13 => false,
    :R_se => true,
    :T_e => false
    ]),
    :J1_4)


## Geometric Parameters
L = 0.2
A = 1.0
Δ = 0.1
## Thermo=hysical Properties
λ = 0.963
c = 650
ρ = 1300
## Inside and Outside Conditions
R_si = 0.13
R_se = 0.04
T_i = 21.0
T_i = 294.15
T_e = 0.0
T_e = 273.15
T_0 = 12.0
T_0 = 285.15
generate_model!(wall)
simplify_model!(wall)
## Set Parameters

p = [
    wall.elements[:T_i].sys.Se => T_i,
    wall.elements[:T_e].sys.Se => T_e,
    wall.elements[:R_si].sys.R => R_si,
    wall.elements[:R_se].sys.R => R_se,
    wall.elements[:R_l].sys.R => Δ/λ/A,
    wall.elements[:R_ll].sys.R => Δ/λ/A,
    wall.elements[:C_l].sys.C => ρ*Δ*A*c,
    wall.elements[:C_ll].sys.C => ρ*Δ*A*c,
    wall.elements[:C_lll].sys.C => ρ*Δ*A*c
    ]|>Dict
ps = ArrayPartition(map(x -> [p[x]], parameters(wall.model))...)

u0 = [
    ρ*Δ*A*c*T_0,
    ρ*Δ*A*c*T_0,
    ρ*Δ*A*c*T_0
    ]
u0 = ArrayPartition(map(x -> [u0[x]], states(wall.model))...)

tspan = (0.0, 40.0*3600.0)
func = build_torn_function(wall.model)
prob = ODEProblem(func, u0, tspan, ps)
sol = solve(prob, Tsit5(), reltol = 1e-5)
plot(sol.t, sol[wall.elements[:R_l].sys.e])

## Geometric Parameters
L = 0.2
A = 1.0
Δ = 0.1
## Thermo=hysical Properties
λ = 0.963
c = 650
ρ = 1300
## Inside and Outside Conditions
R_si = 0.13
R_se = 0.04
T_i = 294.15
T_e = 273.15
T_0 = 285.15
R_l = Δ/λ/A
C_l = ρ*Δ*A*c
u0 = [
    ρ*Δ*A*c*T_0,
    ρ*Δ*A*c*T_0,
    ρ*Δ*A*c*T_0
    ]
function sys!(du, u, t, p)
    du[1] = -(1/R_si+1/(R_l+C_l))*u[1]+1/(R_l*C_l)*u[1]+1/R_si*T_i
    du[2] = 1/(R_l*C_l)*u[1]-(1/R_l/C_l+1/R_l/C_l)*u[2]+1/R_l/C_l*u[3]
    du[3] = 1/R_l/C_l/u[2]-(1/R_l/C_l+1/R_se/C_l)*u[3]-1/R_se*T_e
end

prob = ODEProblem(sys!, u0, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-5)
TWi = getindex.(sol.u,1)./C_l
TW_mid = getindex.(sol.u,2)./C_l
TWe = getindex.(sol.u,3)./C_l
ΔT = T_0-T_e
plot(sol.t, TWi .- T_e/ΔT)