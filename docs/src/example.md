# Package Imports

## Mass - Spring - Damper Example

```@example 1
using SystemBondGraphs
using ModelingToolkit: inputs
using OrdinaryDiffEq
using CairoMakie
using Latexify

set_theme!(theme_latexfonts())
```

Create the time variable and the empty bond graph:

```@example 1
@variables t
bg = BondGraph(t)
```

Elements in the bond graph are added with their appropriate add function and assigned a unique name in the bond graph:

```@example 1
add_R!(bg, :R)
add_C!(bg, :C)
add_I!(bg, :I)
add_Se!(bg, :Se)
```

Next, the 1-Junction is added:

```@example 1
add_1J!(bg, :J1_1)
```

The edges are connected with the `add_bond!` function which modifies the bond graph to have a power bond connecting the nodes. The function is takes the source and sink edge for the power bond. 

```@example 1
add_bond!(bg, :J1_1, :R, :e1)
add_bond!(bg, :J1_1, :C, :e2)
add_bond!(bg, :J1_1, :I, :e3)
add_bond!(bg, :Se, :J1_1, :e4)
```

The system of equations is equation from the bond graph with the `generate_model` function and then simplified with structural_simplify from `ModelingToolkit.jl`:

```@example 1
sys = generate_model(bg)
sys, _ = structural_simplify(sys, (inputs(sys), []))
```

To simulate the system, the solution requires inital conditions, parameters, and timespan to be set. See DifferentialEquations.jl for more information on the methods for solving the ODESystem.

```@example 1
u0 = [
    bg[:C].model.q => 10.0,
    bg[:I].model.p => -1.0
]
ps = [
    bg[:R].model.R => 1.0,
    bg[:C].model.C => 0.01,
    bg[:I].model.I => 1.0,
    bg[:Se].model.Se => 0.0
    ]
tspan = (0.0, 10.0)
```

Finally, solving and plotting the solution gives:

```@example 1
prob = ODEProblem(sys, u0, tspan, ps)
sol = solve(prob, Tsit5())
f, ax, p = CairoMakie.plot(sol,  axis=(xlabel="Time (t)", ))
save("msd.png", f) # hide
```

![Mass Spring Damper Image](msd.png)
