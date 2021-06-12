using Pkg
Pkg.activate("venv")
using BondGraphs
using ModelingToolkit
using DifferentialEquations
using Plots
using Symbolics
using SymbolicUtils
using Unitful
using MonteCarloMeasurements
using Turing
using UnitfulRecipes
using StatsPlots
display("Imports Finished")
## ## ## ## ## ## ## ## ## ## ## ## ## ## 
##  Linear Mass Spring Damper Systems  ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## 
display("Defining System")
@variables t α
msd = BondGraph(t) 
# Elements
add_R!(msd, :r1) 
add_C!(msd, :c1)
add_I!(msd, :i1)
F(t, p) = sin(p[1] * t) * p[2]
@register F(t, p)
@parameters α ω ωc
add_Se!(msd, F, [ω, α], :se)
# Junctions
add_1J!(msd, Dict([
    :r1 => true, 
    :c1 => true, 
    :i1 => true, 
    :se => false]),
    :J1)
generate_model!(msd)
simplify_model!(msd)

# u0 = Dict(
#     msd.elements[:i1].sys.p => 0.0u"kg*m/s",
#     msd.elements[:c1].sys.q => 0.0u"m"
#     )
# u0 = ArrayPartition(map(x->[u0[x]], states(msd.model))...)
# ps = Dict(
#     msd.elements[:r1].sys.R => 0.1u"N*s/m",
#     msd.elements[:c1].sys.C=>0.1u"m/N",
#     msd.elements[:i1].sys.I=>1.0u"kg",
#     msd.elements[:se].sys.α => 1.0u"N",
#     msd.elements[:se].sys.ω=>2.0u"rad/s",
#     )
# ps = ArrayPartition(map(x->[ps[x]], parameters(msd.model))...)

# tspan = (0.0u"s", 10.0u"s")

u0 = Dict(
    msd.elements[:i1].sys.p => 0.0,
    msd.elements[:c1].sys.q => 0.0
    )
u0 = ArrayPartition(map(x -> [u0[x]], states(msd.model))...)
ps = Dict(
    msd.elements[:r1].sys.R => 0.1,
    msd.elements[:c1].sys.C => 0.1,
    msd.elements[:i1].sys.I => 1.0,
    msd.elements[:se].sys.α => 1.0 + 0.2Particles(2000),
    msd.elements[:se].sys.ω => 2.0,
    )
ps = ArrayPartition(map(x -> [ps[x]], parameters(msd.model))...)

tspan = (0.0, 10.0)

display("Generating Model")
func = ODEFunction(msd.model)
prob = ODEProblem(func, u0, tspan, ps)
display("Solving")
sol = solve(prob, Tsit5())
display("Plotting")
plot(sol.t, getindex.(sol.u, 2))
mcplot(sol.t, getindex.(sol.u, 2))
errorbarplot(sol.t, getindex.(sol.u, 2))
ribbonplot(sol.t, getindex.(sol.u, 2))