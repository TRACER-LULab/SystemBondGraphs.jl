using SystemBondGraphs
using ModelingToolkit: linearize_symbolic, inputs
# Create Bond Graphs and Model
@variables t ω(t)
mat = BondGraph(t)
# Add One-Ports
add_I!(mat, :L)
add_R!(mat, :R)
add_I!(mat, :J)
add_R!(mat, :kf)
# Add Inputs
add_Se!(mat, :τ)
add_Se!(mat, :Va)
# Add Multiports
add_GY!(mat, :km)
# Add Junctions
add_1J!(mat, :J11)
add_1J!(mat, :J12)
# Add bonds
add_bond!(mat, :Va, :J11, :edge_1)
add_bond!(mat, :J11, :L, :edge_2)
add_bond!(mat, :J11, :R, :edge_3)
add_bond!(mat, :J11, :km, :edge_4)
add_bond!(mat, :km, :J12, :edge_5)
add_bond!(mat, :J12, :J, :edge_6)
add_bond!(mat, :J12, :kf, :edge_7)
add_bond!(mat, :τ, :J12, :edge_8)
# Generate Model
model = generate_model(mat)
parameters(model)
## Generate the State-Space Model
(; A, B, C, D), sys =
    linearize_symbolic(model, inputs(model), [mat[:J].model.f], simplify = true)
A * unknowns(sys) + B * inputs(sys)
C * unknowns(sys) + D * inputs(sys)
