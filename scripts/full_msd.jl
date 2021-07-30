using DrWatson
@quickactivate "BondGraphModeling"
using Pkg
Pkg.activate(".")
Pkg.instantiate()
##
using BondGraphs
using ModelingToolkit
using DifferentialEquations
using Plots
using Symbolics
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
@parameters Se
add_Se!(msd, Se, :se)
# Junctions
add_1J!(msd, Dict([
    :r1 => true, 
    :c1 => true, 
    :i1 => true, 
    :se => false]),
    :J1)
generate_model!(msd)
simplify_model!(msd)

ps = Dict(
    msd.elements[:r1].sys.R => 1.0,
    msd.elements[:c1].sys.C => 0.1,
    msd.elements[:i1].sys.I => 1.0
    )
    
C = [1 1; 1 1]
D = [0; 0]

ẋ = Ax + Bu
y = Cx + Du
tf = transfer_function(msd, ps, C, D);
G = eval(tf[1])
abs.(G(im))
ω = 10.0.^range(-4, 6.0, length=10000)
plot(ω, transpose(abs.(reduce(hcat, (G.(ω * im))))), xscale=:log10, yscale=:log10)
##
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
# mcplot(sol.t, getindex.(sol.u, 2))
# errorbarplot(sol.t, getindex.(sol.u, 2))
# ribbonplot(sol.t, getindex.(sol.u, 2))

## Transverse HalfCar TFs
@variables t e(t) f(t) 
halfcar = BondGraph(t)
@variables Sf1
add_Sf!(halfcar, Sf1, :sf1)
add_C!(halfcar, :c2)
add_Bond!(halfcar, :b3)
add_0J!(halfcar, Dict([
    :sf1 => false,
    :c2 => true,
    :b3 => true
    ]),
    :J1)
add_I!(halfcar, :i4)
add_Se!(halfcar, :se5)
add_Bond!(halfcar, :b6)
add_Bond!(halfcar, :b24)
add_Bond!(halfcar, :b23)
add_1J!(halfcar, Dict([
    :b3 => false,
    :i4 => true,
    :se5 => false,
    :b6 => true,
    :b23 => true,
    :b24 => true
    ]), 
    :J2)
add_Bond!(halfcar, :b7)
add_Bond!(halfcar, :b10)
add_0J!(halfcar, Dict([
    :b6 => false,
    :b7 => true,
    :b10 => true
    ]),
    :J3)
add_C!(halfcar, :c8)
add_R!(halfcar, :r9)
add_1J!(halfcar, Dict([
    :b7 => false,
    :c8 => true,
    :r9 => true
    ]),
    :J4)
add_Se!(halfcar, :se11)
add_I!(halfcar, :i12)
add_Bond!(halfcar, :b13)
add_1J!(halfcar, Dict([
    :b10 => false,
    :se11 => true,
    :i12 => true,
    :b13 => false
    ]),
    :J5)
add_Bond!(halfcar, :b14)
add_Bond!(halfcar, :b17)
add_0J!(halfcar, Dict([
    :b13 => true,
    :b14 => true,
    :b17 => false
    ]),
    :J6)
add_C!(halfcar, :c15)
add_R!(halfcar, :r16)
add_1J!(halfcar, Dict([
    :b14 => false,
:c15 => true,
    :r16 => true
    ]),
    :J7)
add_Se!(halfcar, :se18)
add_I!(halfcar, :i19)
add_C!(halfcar, :c20)
add_Bond!(halfcar, :b21)
add_Bond!(halfcar, :b27)
add_1J!(halfcar, Dict([
    :b17 => true,
    :se18 => false,
    :i19 => true,
    :c20 => true,
    :b21 => false,
    :b27 => false
    ]),
    :J8)
add_C!(halfcar, :c22)
add_0J!(halfcar, Dict([
    :b21 => true,
    :c22 => true,
    :b23 => false
    ]),
    :J9)
add_R!(halfcar, :r25)
add_0J!(halfcar, Dict([
    :b24 => false,
    :r25 => true,
    :b27 => true
    ]),
    :J10)