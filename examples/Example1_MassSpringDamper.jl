using BondGraphs

# Create a bond graph
@variables t
bg = BondGraph(t, :msd)

# Create the Elements in the bond graph
add_R!(bg, :R)
add_C!(bg, :C)
add_I!(bg, :I)
add_Se!(bg, :Se)

# Add the Junction
add_1J!(bg, :J1_1)

# Connect the nodes
add_bond!(bg, :R, :J1_1, :edge_1)
add_bond!(bg, :C, :J1_1, :edge_2)
add_bond!(bg, :J1_1, :I, :edge_3)
add_bond!(bg, :Se, :J1_1, :edge_4)

# Form the sytem of equations
sys = generate_model(bg)
sys = structural_simplify(sys)
