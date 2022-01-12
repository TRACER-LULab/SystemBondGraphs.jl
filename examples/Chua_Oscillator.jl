using BondGraphs
using ModelingToolkit
using DifferentialEquations
## --- 
@variables t
chua = BondGraph(t)
# Add Nonlinear element
@parameters Ga Gb E
Φr(e, f, t, p) = p[2] * e + 1 / 2 * (p[2] - p[1]) * (abs(e + p[3]) - abs(e - p[3]))
add_R!(chua, :f => Φr, [Ga, Gb, E], :Gn)
# Add Linear Elements
add_C!(chua, :C1)
add_C!(chua, :C2)
add_R!(chua, :R)
add_R!(chua, :G)
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
add_1J!(chua, Dict(
  :b8 => true, 
  :C1 => false,
  :Gn => false
), :J14)
## Generate Model
sys = generate_model(chua)
sys = structural_simplify(sys)
# states(sys)
# # Set Parameters
p = [
  chua[:R].R => 1.75,
  chua[:C1].C => 10e-9,
  chua[:C2].C => 100e-9,
  chua[:L].I => 18e-3,
  chua[:R_nonlinear].E => 1.0,
  chua[:R_nonlinear].Ga => -757e-6,
  chua[:R_nonlinear].Gb => -409e-6
] |> Dict
u0 = [
  chua[:C1].q => -1.35164*p[chua[:C1].C],
  chua[:C2].q => 9.13959*p[chua[:C2].C],
  chua[:L].p  => 9.2869*p[chua[:L].I],
] |> Dict
tspan = (0.0, 0.05)
prob = ODAEProblem(sys, u0, tspan, p)
sol = solve(prob, Rosenbrock23())