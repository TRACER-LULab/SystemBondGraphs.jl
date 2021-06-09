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
r(t, p) = tanh(t)
@register r(t, p)
add_Se!(msd, r, [], :se)
# Junctions
add_1J!(msd, Dict([
    :r1 => true, 
    :c1 => true, 
    :i1 => true, 
    :se => false]),
    :J1)
u0 = [
    msd.elements[:c1].sys.q => 0.0u"m",
    msd.elements[:i1].sys.p => 0.0u"N*m/s",
    ]
ps = [
    # msd.elements[:se].sys.α => 1.0u"N",
    msd.elements[:c1].sys.C => 0.1u"m/N"
    msd.elements[:r1].sys.R => 10.0u"N*s/m",
    msd.elements[:i1].sys.I => 1.0u"kg"
]
# Set TimeSpan
tspan = (0.0u"s", 10.0u"s")
display("Generating Model")
generate_model!(msd)
simplify_model!(msd)
prob = ODEProblem(msd.model, u0, tspan, ps)
display("Solving")
sol = solve(prob, Tsit5())
display("Plotting")
plot(sol)

function sys!(du, u, p, t)
    c, r, i = p
    du[1] = tanh(t) - u[2] / c - r * u[1] / i
    du[2] = u[1] / i
end

u0 = [
    0.0u"N/s",
    0.0u"m"
]
ps = [
    0.1u"m/N",
    10.0u"N*s/m",
    1.0u"kg"
]
prob = ODEProblem(sys!, u0, (0.0u"s", 10.0u"s"), ps)
solve(prob)

g = (y, p, t) -> 0.5 * y
u0 = 1.5±0.2
prob = ODEProblem(g, u0, (0.0, 1.0))
sol = solve(prob, Tsit5())

display(plot(sol))