#=
Bond Graph System from:
> Medina H, Carpenter N. Transverse Third-Damper to Improve Suspension Systems for In-Wheel Motor Driven Electric Vehicles. Journal of Dynamic Systems, Measurement, and Control. 2023 Mar 1;145(3):031003.
=#
using BondGraphs
using DifferentialEquations
using Plots
using ControlSystemsMTK, ControlSystems
## Analysis for Third Damper Model
@parameters t
# Create Empty Bondgraph
conventional = BondGraph(t, :car)
# Add Flow Sources
add_Sf!(conventional, :Vin)
# Add Effort Sources
add_Se!(conventional, :m_us1_g)
add_Se!(conventional, :m_s_g)
add_Se!(conventional, :m_us2_g)
# Add Compliance Ports
add_C!(conventional, :kt1)
add_C!(conventional, :k1)
add_C!(conventional, :k2)
add_C!(conventional, :kt2)
add_C!(conventional, :k)
# Add Resistive Ports
add_R!(conventional, :b1)
add_R!(conventional, :b2)
# Add Inertial Ports
add_I!(conventional, :m_us1)
add_I!(conventional, :m_s)
add_I!(conventional, :m_us2)
# Add 1 Junctions
add_1J!(conventional, :J1_1)
add_1J!(conventional, :J1_2)
add_1J!(conventional, :J1_3)
add_1J!(conventional, :J1_4)
add_1J!(conventional, :J1_5)
# Add 0 Junctions
add_0J!(conventional, :J0_1)
add_0J!(conventional, :J0_2)
add_0J!(conventional, :J0_3)
add_0J!(conventional, :J0_4)
# Add Bonds
add_bond!(conventional, :Vin, :J0_1, :e1)
add_bond!(conventional, :J0_1, :kt1, :e2)
add_bond!(conventional, :J0_1, :J1_1, :e3)
add_bond!(conventional, :J1_1, :m_us1, :e4)
add_bond!(conventional, :J1_1, :m_us1_g, :e5)
add_bond!(conventional, :J1_1, :J0_2, :e6)
add_bond!(conventional, :J1_1, :J0_4, :e23)
add_bond!(conventional, :J0_2, :J1_2, :e7)
add_bond!(conventional, :J1_2, :k1, :e8)
add_bond!(conventional, :J1_2, :b1, :e9)
add_bond!(conventional, :J0_2, :J1_3, :e10)
add_bond!(conventional, :J1_3, :m_s_g, :e11)
add_bond!(conventional, :J1_3, :m_s, :e12)
add_bond!(conventional, :J0_3, :J1_3, :e13)
add_bond!(conventional, :J0_3, :J1_4, :e14)
add_bond!(conventional, :J1_4, :k2, :e15)
add_bond!(conventional, :J1_4, :b2, :e16)
add_bond!(conventional, :J1_5, :J0_3, :e17)
add_bond!(conventional, :J1_5, :m_us2_g, :e18)
add_bond!(conventional, :J1_5, :m_us2, :e19)
add_bond!(conventional, :J1_5, :kt2, :e20)
add_bond!(conventional, :J0_4, :J1_5, :e21)
add_bond!(conventional, :J0_4, :k, :e22)
# TRANSVERSE-DAMPER Model
third_damper = copy(conventional)
# Add Transverse Damper
add_R!(third_damper, :b)
add_0J!(third_damper, :J0_5)
add_bond!(third_damper, :J1_1, :J0_5, :e24)
add_bond!(third_damper, :J0_5, :J1_5, :e27)
# PLOTTING
p = plot(layout=(2,1), size = (800,800))
setPlotScale("dB")
for m in [40, 60, 80]
    # Generate Models
    sys = generate_model(conventional)
    # Parameters
    b = 1964.0
    subs_dict = [
        conventional[:m_us1_g].Se => conventional[:m_us1].I * 9.81,
        conventional[:m_us2_g].Se => conventional[:m_us2].I * 9.81,
        conventional[:m_s_g].Se => conventional[:m_s].I * 9.81
    ] |> Dict
    sys = substitute(sys, subs_dict)
    #
    op = [
        conventional[:m_s].I => 680,
        conventional[:m_us1].I => m,
        conventional[:m_us2].I => m,
        conventional[:k1].C => 1 / 31.5e3,
        conventional[:k2].C => 1 / 31.5e3,
        conventional[:kt1].C => 1 / 220e3,
        conventional[:kt2].C => 1 / 220e3,
        conventional[:k].C => 1 / 220e3,
        conventional[:b1].R => b,
        conventional[:b2].R => b,
        conventional[:Vin].Sf => 0.0,
        conventional[:kt1].q => 0.0,
        conventional[:k1].q => 0.0,
        conventional[:k2].q => 0.0,
        conventional[:kt2].q => 0.0,
        conventional[:k].q => 0.0,
        conventional[:m_us1].p => 0.0,
        conventional[:m_s].p => 0.0,
        conventional[:m_us2].p => 0.0
    ] |> Dict
    state_space = named_ss(sys, [conventional[:Vin].Sf], [conventional[:m_us2].e], simplify=true, op=op)
    bodeplot!(state_space, hz=true, scale=:identity, xlims=[0, 30], xticks=0:5:30, label="Conventional: Mass = $m")
    # -------------------------------------------------
    # -------------------------------------------------
    # -------------------------------------------------

    add_bond!(third_damper, :J0_5, :b, :e25)
    new_sys = generate_model(third_damper)
    new_p = [
        third_damper[:m_us1_g].Se => third_damper[:m_us1].I * 9.81,
        third_damper[:m_us2_g].Se => third_damper[:m_us2].I * 9.81,
        third_damper[:m_s_g].Se => third_damper[:m_s].I * 9.81
    ] |> Dict
    new_sys = substitute(new_sys, new_p)
    new_op = [
        third_damper[:m_s].I => 680,
        third_damper[:m_us1].I => m,
        third_damper[:m_us2].I => m,
        third_damper[:k1].C => 1 / 31.5e3,
        third_damper[:k2].C => 1 / 31.5e3,
        third_damper[:kt1].C => 1 / 220e3,
        third_damper[:kt2].C => 1 / 220e3,
        third_damper[:k].C => 1 / 220e3,
        third_damper[:b1].R => 1964.0,
        third_damper[:b2].R => 1964.0,
        third_damper[:Vin].Sf => 0.0,
        third_damper[:kt1].q => 0.0,
        third_damper[:k1].q => 0.0,
        third_damper[:k2].q => 0.0,
        third_damper[:kt2].q => 0.0,
        third_damper[:k].q => 0.0,
        third_damper[:m_us1].p => 0.0,
        third_damper[:m_s].p => 0.0,
        third_damper[:m_us2].p => 0.0,
        third_damper[:b].R => 1964
    ] |> Dict
    new_state_space = named_ss(new_sys, [third_damper[:Vin].Sf], [third_damper[:m_us2].e], simplify=true, op=new_op)
    bodeplot!(new_state_space,  hz=true, scale=:identity, xlims=[0, 30], xticks=0:5:30, label="Third Damper: Mass = $m", linestyle=:dash, )
end
p
