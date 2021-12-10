## Package Imports
using BondGraphs
using DifferentialEquations
using ModelingToolkit
using Plots
## Create BondGraph Object
@parameters t
sys = BondGraph(t)
## Add Elements
# Input Velocity
add_Sf!(sys, :Sf)
# Control Velocity
add_Se!(sys, :Se)
# Add Masses
add_I!(sys, :m)
add_I!(sys, :mc)
# Add Dampers
add_R!(sys, :b)
add_R!(sys, :bc)
# Add Springs
add_C!(sys, :k)
add_C!(sys, :kc)
# Add Resistor
add_R!(sys, :R)
# Add Miscellaneous Bonds
for b in Symbol.("b".*["2","5","7","9","12","13"])
    add_Bond!(sys, b)
end
## Add Two-Ports
@parameters τc
add_GY!(sys, τc, :b13, :b12, :τc)
## Add Multi-Ports
add_1J!(sys, Dict([
    :b => false,
    :k => false,
    :b2 => true
    ]), :J11)
add_0J!(sys, Dict([
    :Sf => true,
    :b2 => false,
    :b5 => false
    ]), :J01)
add_1J!(sys, Dict([
    :b5 => true,
    :m => false,
    :b7 => false
    ]), :J12)
add_0J!(sys, Dict([
    :b7 => true,
    :b9 => false,
    :mc => false
    ]), :J02)
add_1J!(sys, Dict([
    :b9 => true,
    :bc => false,
    :b12 => true,
    :kc => false
    ]), :J13)
add_1J!(sys, Dict([
    :Se => true,
    :R => false,
    :b13 => false
    ]), :J14)
## Create the model
model = generate_model(sys)
model = structural_simplify(model)
## Generate Transfer Functions
generate_model!(sys)
sys.model = structural_simplify(sys.model)
state_space(sys)