# # Nonlinear Chaotic Chua Oscillator
using BondGraphs
using ModelingToolkit
using DifferentialEquations
# # Make the BondGraph
@variables t
chua = BondGraph(t, :chua)
# # Add Nonlinear Resistor Equation
@parameters Ga Gb E
Φr(e, f, t, p) = p[2] * e + 1 / 2 * (p[1] - p[2]) * (abs(e + p[3]) - abs(e - p[3]))
add_R!(chua, :f => Φr, [Ga, Gb, E], :Gn)
# # Add Linear Elements
add_R!(chua, :G)
add_R!(chua, :R)
add_C!(chua, :C1)
add_C!(chua, :C2)
add_I!(chua, :L)
# # Add Bonds
add_Bond!(chua, :b2)
add_Bond!(chua, :b3)
add_Bond!(chua, :b5)
add_Bond!(chua, :b8)
# # Add Junctions
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
        :L => false,
        :R => false,
        :b8 => false
    ), :J13)
add_0J!(chua, Dict(
        :b8 => true,
        :C1 => false,
        :Gn => false
    ), :J02)
# # Generate the model and solve
sys = generate_model(chua)
model = structural_simplify(sys, simplify=true)
u0 = [
    chua[:C1].q => 1 / 10 / 9.13959,
    chua[:C2].q => 2 / -1.35164,
    chua[:L].p => 1 / 7 / -59.2869
]
ps = [
    chua[:R].R => 1 / 0.7,
    chua[:G].R => 1 / 0.7,
    chua[:C1].C => 1 / 10.0,
    chua[:C2].C => 1 / 0.5,
    chua[:L].I => 1 / 7.0,
    chua[:Gn].E => 1.0,
    chua[:Gn].Ga => 10.0,
    chua[:Gn].Gb => 100.0
]
tspan = (0.0, 10000.0)
prob = ODEProblem(model, u0, tspan, ps)
sol = solve(prob, Tsit5())
# # Plot the Results
plot(sol, vars = (chua[:C1].q, chua[:C2].q, chua[:L].p))
