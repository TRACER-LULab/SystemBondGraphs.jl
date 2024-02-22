## Imports
using BondGraphs
using DifferentialEquations
# Create BondGraph
@variables t
bg = BondGraph(t, :bg)
# Add One Ports
add_Se!(bg, :in)
add_I!(bg, :mc)
add_I!(bg, :mpx)
add_I!(bg, :J)
add_Se!(bg, :mpg)
add_I!(bg, :mpy)
# Add Transformers
@variables θ(t) x(t)
@parameters l
add_MTF!(bg, cos(θ) * l, :v_x)
add_MTF!(bg, -sin(θ) * l, :v_y)
# Add Junctions
add_1J!(bg, :v_c_x)
add_0J!(bg, :J0_x)
add_1J!(bg, :v_p_x)
add_1J!(bg, :v_p_ω)
add_0J!(bg, :J0_y)
add_1J!(bg, :v_p_y)
# Add Edges
add_bond!(bg, :v_p_ω, :v_x, :edge_1)
add_bond!(bg, :v_x, :J0_x, :edge_2)
add_bond!(bg, :v_p_ω, :v_y, :edge_3)
add_bond!(bg, :v_y, :J0_y, :edge_4)
add_bond!(bg, :in, :v_c_x, :edge_5)
add_bond!(bg, :v_c_x, :mc, :edge_6)
add_bond!(bg, :J0_x, :v_c_x, :edge_7)
add_bond!(bg, :v_p_x, :J0_x, :edge_8)
add_bond!(bg, :mpx, :v_p_x, :edge_9)
add_bond!(bg, :v_p_ω, :J, :edge_10)
add_bond!(bg, :J0_y, :v_p_y, :edge_11)
add_bond!(bg, :mpg, :v_p_y, :edge_12)
add_bond!(bg, :v_p_y, :mpy, :edge_13)
# Add the MTF Equations
D = Differential(bg.model.iv)
eqns = [
    D(θ) ~ bg[:J].f,
    D(x) ~ bg[:mc].f,
]
@named theta_model = ODESystem(eqns, bg.model.iv, [θ, x], [])
# Extend and build the model
sys = generate_model(bg)
sys = extend(sys, theta_model)
sys = structural_simplify(sys)
# Initialize and test the system
u0 = [
    x => 0.0,
    θ => -5 * 3.14159 / 180,
    bg[:mpx].f => 0.0,
    bg[:v_y].f_out => 0.0,
]
p = [
    l => 1.0,
    bg[:mc].I => 1.0,
    bg[:mpx].I => 1.0,
    bg[:mpy].I => 1.0,
    bg[:J].I => 1.0^2 * 1.0,
    bg[:mpg].Se => 1.0 * 9.81,
    bg[:in].Se => 0.0,
]
tspan = (0.0, 10.0)
new_prob = ODAEProblem(sys, u0, tspan, p)
sol = solve(new_prob, Tsit5())
