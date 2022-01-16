using BondGraphs
using ModelingToolkit
using DifferentialEquations
## --- 
@variables t
chua = BondGraph(t)
# Add Nonlinear element
@parameters Ga Gb E
Φr(e, f, t, p) = p[2] * e + 1 / 2 * (p[1] - p[2]) * (abs(e + p[3]) - abs(e - p[3]))
add_R!(chua, :f => Φr, [Ga, Gb, E], :Gn)
# Add Linear Elements
add_R!(chua, :G)
add_R!(chua, :R)
add_C!(chua, :C1)
add_C!(chua, :C2)
add_I!(chua, :L)
# Add Bonds
add_Bond!(chua, :b2)
add_Bond!(chua, :b3)
add_Bond!(chua, :b5)
add_Bond!(chua, :b8)
# Add Junctions
add_1J!(chua, Dict(
  :C2 => false,
  :b2 => false
  ), :J11)
add_1J!(chua, Dict(
  :G => false,
  :b3 => false
  ), :J12)
add_0J!(chua, Dict(
  :b2 => true,
  :b3 => true,
  :b5 => false
  ), :J01)
add_1J!(chua, Dict(
  :b5 => true,
  :L  => false,
  :R  => false,
  :b8 => false
  ), :J13)
add_0J!(chua, Dict(
  :b8 => true,
  :C1 => false,
  :Gn => false
  ), :J02)
sys = generate_model(chua)
sys = structural_simplify(sys, simplify = true)
# equations(sys)
# ## # Set Parameters
# p = [
#   chua[:R].R => 1 / 0.7,
#   chua[:C1].C => 1 / 10.0,
#   chua[:C2].C => 1 / 0.5,
#   chua[:L].I => 1 / 7.0,
#   chua[:R_nonlinear].E => 1.0,
#   chua[:R_nonlinear].Ga => 1/4.0,
#   chua[:R_nonlinear].Gb => 1/0.1
# ] |> Dict
# u0 = [
#   chua[:C1].q => -1.35164 * p[chua[:C1].C],
#   chua[:C2].q => 9.13959 * p[chua[:C2].C],
#   chua[:L].p => 9.2869 * p[chua[:L].I],
# ] |> Dict
# tspan = (0.0, 10.0)
# prob = ODAEProblem(sys, u0, tspan, p)
# sol = solve(prob, Rosenbrock23())
# ##
# plot(sol[chua[:C1].e], sol[chua[:C2].e])
# ##
# # ---- traditional implementation of chua's circuit ----
# @variables t v1(t) v2(t) i3(t)
# @parameters c1 c2 L Ga Gb E G
# D = Differential(t)
# f = IfElse.ifelse(-1<v1<1, -Gb*v1, IfElse.ifelse(v>1, -Gb-Ga*(v1-1.0), Gb-Ga*(v+1.0)))
# eqns = [
#   D(v1) ~ 1 / c1 * (G * (v2 - v1)) - ((-1<v1<1)*(-Gb*v)+(v1 >1)*(-Gb-Ga*(v1-1.0))+(v1<-1.0)*(Gb-Ga*(v1+1.0)))
#   D(v2) ~ 1 / c2 * (G * (v1 - v2) + i3)
#   D(i3) ~ -v2 / L
# ]
# # -------------------------
# @named sys = ODESystem(eqns, t)
# # -------------------------
# u0 = [
#   v1 => 9.13959, 
#   v2 => -1.35164, 
#   i3 => -59.2869
# ]

# p = [
#   G => 0.7,
#   c1 => 1 / 10.0,
#   c2 => 1/ 0.5,
#   L => 1 / 7.0,
#   E => 1.0,
#   Ga => 10.0,
#   Gb => 100.0
# ]
# # --------------------------
# prob = ODEProblem(sys, u0, (0, 600.0), p)
# sol = solve(prob, RK4(), dt = 0.02)
# plot(sol, vars = [(v1,v2,i3)])
# plot(sol, vars = [(v1, v2)])
# plot(sol, vars = [(v1, i3)])
# plot(sol, vars = [(v2, i3)])