# Package Imports

## Mass - Spring - Damper Example

```@example 1
using OrdinaryDiffEq
using SystemBondGraphs
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

For fiting multiple models to the same dataset, 

```@example 1
models = Dict(
    Gent => ComponentVector(μ=240e-3, J_m=80.0),
    EdwardVilgis => ComponentVector(Ns=0.10, Nc=0.20, α=0.001, η=0.001),
    NeoHookean => ComponentVector(μ=200e-3),
    Beda => ComponentVector(C1=0.1237, C2=0.0424, C3=7.84e-5, K1=0.0168, α=0.9, β=0.68, ζ=3.015)
)

sol = Dict{Any, SciMLSolution}()
for (ψ, p_0) in models
    HEProblem = HyperelasticProblem(ψ(), treloar_data, p_0,  ad_type = AutoForwardDiff())
    sol[ψ] = solve(HEProblem, NelderMead())
end
return sol # hide
```

To predict the reponse of a model to a proivded dataset and parameters, a `predict` function is provided:

```@example 1
f = Figure(size = (800,800))
ax = Makie.Axis(f[1,1], xlabel = "Stretch [-]", ylabel = "Stress [kg/cm²]")
for (ψ, p) in sol
    pred = predict(ψ(), treloar_data, p.u, ad_type = AutoForwardDiff())
    lines!(ax, getindex.(pred.data.λ, 1), getindex.(pred.data.s, 1), label=string(ψ))
end
scatter!(ax, getindex.(treloar_data.data.λ, 1), getindex.(treloar_data.data.s, 1), label = "Treloar 1944 Experimental", color = :black)
axislegend(position = :lt)
save("treloar_data_fits.png", f) # hide
```

![](treloar_data_fits.png)

While the majority of the models provided by `Hyperelastics.jl` are based on closed form strain energy density functions, a selection of data-driven models are proivded. For example, the `SussmanBathe` model is created with:

```@example 1
using DataInterpolations
ψ = SussmanBathe(treloar_data, k=4, interpolant = QuadraticSpline)
λ₁ = range(extrema(getindex.(treloar_data.data.λ, 1))..., length = 100)
uniaxial_prediction = HyperelasticUniaxialTest(λ₁, name = "Prediction")
pred = predict(ψ, uniaxial_prediction, [])
λ₁ = getindex.(treloar_data.data.λ, 1)
s₁ = getindex.(treloar_data.data.s, 1)
λ̂₁ = getindex.(pred.data.λ, 1)
ŝ₁ = getindex.(pred.data.s, 1)

f = Figure(size = (800,800))
ax = Makie.Axis(f[1,1], xlabel = "Stretch [-]", ylabel = "Stress [kg/cm²]")
lines!(
    ax, 
    λ̂₁, 
    ŝ₁, 
    label = "Sussman-Bathe Approximation"
)

scatter!(
        ax,
        λ₁, 
        s₁, 
        label = "Treloar 1944 Experimental",
        color = :black
    )
axislegend(position = :lt)
save("sussman_bathe.png", f) # hide
```
![](sussman_bathe.png)