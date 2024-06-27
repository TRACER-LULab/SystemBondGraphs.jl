using SystemBondGraphs
using ModelingToolkit: linearize_symbolic, inputs
using ControlSystems
# Create Bond Graphs and Model
@variables t ω(t)
bg = BondGraph(t)
# Add One-Ports
add_I!(bg, :L)
add_R!(bg, :R)
add_I!(bg, :J)
add_R!(bg, :kf)
# Add Inputs
add_Se!(bg, :τ)
add_Se!(bg, :Va)
# Add Multiports
add_GY!(bg, :km)
# Add Junctions
add_1J!(bg, :J11)
add_1J!(bg, :J12)
# Add bonds
add_bond!(bg, :Va, :J11, :edge_1)
add_bond!(bg, :J11, :L, :edge_2)
add_bond!(bg, :J11, :R, :edge_3)
add_bond!(bg, :J11, :km, :edge_4)
add_bond!(bg, :km, :J12, :edge_5)
add_bond!(bg, :J12, :J, :edge_6)
add_bond!(bg, :J12, :kf, :edge_7)
add_bond!(bg, :τ, :J12, :edge_8)
# Generate Model
sys = generate_model(bg)
sys_ss, _ = structural_simplify(sys, (inputs(sys), []))
parameters(sys)
## Generate the State-Space Model
(; A, B, C, D), sys =
    linearize_symbolic(sys, inputs(sys), [bg[:J].model.f], simplify=true)
A * unknowns(sys) + B * inputs(sys)
C * unknowns(sys) + D * inputs(sys)
## Generate Transfer Function
ps = [
    bg[:km].model.r => 1.0,
    bg[:L].model.I => 1.0,
    bg[:R].model.R => 1.0,
    bg[:J].model.I => 1.0,
    bg[:kf].model.R => 1.0,
    bg[:τ].model.Se => 0.0,
    bg[:Va].model.Se => 0.0,
]
A, B, C, D = map(x->Symbolics.value.(substitute.(x, (ps,))), [A, B, C, D])
state_space = ss(A, B, C, D)
tf(state_space)
