#=
Bond Graph System from:
> Medina H, Carpenter N. Transverse Third-Damper to Improve Suspension Systems for In-Wheel Motor Driven Electric Vehicles. Journal of Dynamic Systems, Measurement, and Control. 2023 Mar 1;145(3):031003.
=#
using SystemBondGraphs
using OrdinaryDiffEq
using ControlSystems
using ModelingToolkit: inputs, linearize_symbolic
using CairoMakie
# Plot Setup
set_theme!(theme_latexfonts())
f = Figure(size = (1500, 750) ./ 2)
mag_ax = Axis(
    f[1, 1],
    xlabel = "Frequency (Hz)",
    ylabel = L"\text{Amplitude Ratio }\left[\frac{e_{mus_2}}{V_{in} b}\right]",
    xticks = 0:5:30,
    yticks = 0:0.5:3.0,
)
colors = Makie.wong_colors()
linestyle = [:solid, :dash]
mag_plots = []
labels = []

# Create Bondgraph for Conventional Half Car Model
@parameters t
conventional = BondGraph(t)

# Add Flow Sources
add_Sf!(conventional, :Vin)

# Add Effort Sources
add_Se!(conventional, :m_us1_g)
add_Se!(conventional, :m_s_g)
add_Se!(conventional, :m_us2_g)

# Add Compliance Elements
add_C!(conventional, :kt1)
add_C!(conventional, :k1)
add_C!(conventional, :k2)
add_C!(conventional, :kt2)
add_C!(conventional, :k)

# Add Resistive Elements
add_R!(conventional, :b1)
add_R!(conventional, :b2)

# Add Inertial Elements
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

# Create a copy of the tradiational model and add a transverse damper
third_damper = copy(conventional)
# Add Transverse Damper
add_R!(third_damper, :b)
add_0J!(third_damper, :J0_5)
add_bond!(third_damper, :J1_1, :J0_5, :e24)
add_bond!(third_damper, :J0_5, :J1_5, :e27)
add_bond!(third_damper, :J0_5, :b, :e25)

# Set model parameters and substitute
# Parameters
b = 1964.0
ps = Dict([
    conventional[:k1].model.C => 1 / 31.5e3,
    conventional[:k2].model.C => 1 / 31.5e3,
    conventional[:kt1].model.C => 1 / 220e3,
    conventional[:kt2].model.C => 1 / 220e3,
    conventional[:k].model.C => 1 / 220e3,
    conventional[:b1].model.R => b,
    conventional[:b2].model.R => b,
    conventional[:m_s].model.I => 680,
    third_damper[:b].model.R => b,
])

# Operating Point for Linearization
op = Dict([
    conventional[:Vin].model.Sf => 0.0,
    conventional[:kt1].model.q => 0.0,
    conventional[:k1].model.q => 0.0,
    conventional[:k2].model.q => 0.0,
    conventional[:kt2].model.q => 0.0,
    conventional[:k].model.q => 0.0,
    conventional[:m_us1].model.p => 0.0,
    conventional[:m_s].model.p => 0.0,
    conventional[:m_us2].model.p => 0.0,
])
# Generate System of Equations
conventional_model = generate_model(conventional)
conventional_system = substitute(conventional_model, ps)

third_damper_system = substitute(generate_model(third_damper), ps)
systems = Dict("Conventional" => conventional_system, "Third Damper" => third_damper_system)
# Iterate over the different systems
for (i, (name, system)) in enumerate(systems)
    for (j, m) in enumerate([40, 60, 80])
        # Update the mass in the model
        m_ps = Dict([
            conventional[:m_us1_g].model.Se => m * 9.81,
            conventional[:m_us2_g].model.Se => m * 9.81,
            conventional[:m_s_g].model.Se => m * 9.81,
            conventional[:m_us1].model.I => m,
            conventional[:m_us2].model.I => m,
        ])
        system_m = substitute(system, m_ps)
        #
        (; A, B, C, D), _ = linearize_symbolic(
            system_m,
            inputs(system_m),
            [conventional[:m_us2].model.e],
            simplify = true,
            op = op,
        )
        A, B, C, D = map(x -> Symbolics.value.(x), [A, B, C, D])
        state_space = ss(A, B, C, D)
        mag, _, w = bode(state_space, (0.0:0.1:30.0) .* 2π)
        p = lines!(
            mag_ax,
            w ./ 2π,
            mag[1, 1, :] ./ b,
            color = colors[j],
            linestyle = linestyle[i],
        )
        push!(mag_plots, p)
        push!(labels, "$name: $m kg")
    end
end
Legend(f[2, 1], mag_plots, labels, nbanks = 3, tellheight = true, tellwidth = false)
f
