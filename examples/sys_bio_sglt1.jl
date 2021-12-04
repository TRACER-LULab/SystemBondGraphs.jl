using BondGraphs
using ModelingToolkit
using DifferentialEquations
##
@parameters t
sglt = BioBondGraph(t, R=1.0, T=1.0)

# Add Species
add_Ce!(sglt, :Co)
add_Ce!(sglt, :Nao)
add_Ce!(sglt, :CNao)
add_Ce!(sglt, :So)
add_Ce!(sglt, :SCNao)
add_Ce!(sglt, :Ci)
add_Ce!(sglt, :Nai)
add_Ce!(sglt, :CNai)
add_Ce!(sglt, :Si)
add_Ce!(sglt, :SCNai)

# Add Bonds
for i ∈ [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24]
    add_Bond!(sglt, Symbol("B$i"))
end

# Reactions
add_Re!(sglt, :B5, :B6, :r12)
add_Re!(sglt, :B7, :B18, :r25)
add_Re!(sglt, :B10, :B11, :r23)
add_Re!(sglt, :B12, :B13, :r34)
add_Re!(sglt, :B14, :B15, :r45)
add_Re!(sglt, :B19, :B20, :r56)
add_Re!(sglt, :B24, :B1, :r61)

# Add 1 Junctions
add_1J!(sglt, Dict(
    :B2 => true,
    # :B3 => true, 
    :B4 => true, 
    :B5 => false
    ), :J11)
add_1J!(sglt, Dict(
    :B8 => true, 
    :B9 => true, 
    :B10 => false
    ), :J12)
add_1J!(sglt, Dict(
    :B15 => true,
    :B16 => false,
    :B17 => false
    ), :J13)
add_1J!(sglt, Dict(
    :B20 => true, 
    # :B21 => false, 
    :B22 => false, 
    :B23 => false
    ), :J14)

# Add 0 Junctions
add_0J!(sglt, Dict(
    :B1 => true,
    :Co => false, 
    :B2 => false
    ), :j01)
add_0J!(sglt, Dict(
    :Nao => false,
    # :B3 => false, 
    :B4 => false
    ), :j02)
add_0J!(sglt, Dict(
    :CNao => false, 
    :B6 => true,
    :B7 => false,
    :B8 => false
    ), :j03)
add_0J!(sglt, Dict(
    :So => false, 
    :B9 => false
    ), :j04)
add_0J!(sglt, Dict(
    :B11 => true, 
    :B12 => false, 
    :SCNao => false
    ), :j05)
add_0J!(sglt, Dict(
    :SCNai => false,
    :B13 => true, 
    :B14 => false
    ), :j06)
add_0J!(sglt, Dict(
    :B16 => true, 
    :Si => false
    ), :j07)
add_0J!(sglt, Dict(
    :CNai => false,
    :B17 => true, 
    :B18 => true, 
    :B19 => false
    ), :j08)
add_0J!(sglt, Dict(
    :B22 => true, 
    # :B21 => true, 
    :Nai => false
    ), :j09)
add_0J!(sglt, Dict(
    :B23 => true,
    :B24 => false,
    :Ci => false
    ), :j010)
## Create Model
model = generate_model(sglt)
model = structural_simplify(model)

## Log Simplification
using SymbolicUtils
RW = SymbolicUtils.Rewriters
r1 = @acrule log(~x) + log(~y) => log((~x) * (~y))
r2 = @rule log(~x) - log(~y) => log((~x) / (~y))
r3 = @rule (~x) * log(~y) => log((~y)^(~x))
r4 = @rule exp(log(~x)) => ~x
r5 = @acrule exp((~x) + (~y)) => exp(~x) * exp(~y)
rw1 = RW.Fixpoint(RW.Chain([r1, r2, r3, r4, r5]))
rw2 = RW.Prewalk(RW.Chain([r1, r2, r3, r4, r5]))
rw3 = RW.Postwalk(RW.Chain([r1, r2, r3, r4, r5]))
eqns = equations(model)
for i ∈ eachindex(eqns)
    eqns[i] = eqns[i].lhs ~ eqns[i].rhs |> rw3 |> rw2 |> rw1 |> expand
end
defaults = model.defaults
model = ODESystem(eqns, name = :model)
for k ∈ collect(keys(defaults))
    model.defaults[k] = defaults[k]
end
## Set Simulation Conditions
p = [
    sglt[:So].k    => 10.0777,
    sglt[:Si].k    => 10.1247,
    sglt[:Nao].k   => 13.9591,
    sglt[:Nai].k   => 13.9263,
    sglt[:Co].k    => 40.3317,
    sglt[:CNao].k  =>49.1178,
    sglt[:SCNao].k =>0.0989984,
    sglt[:Ci].k    => 0.3457,
    sglt[:CNai].k  => 0.148991,
    sglt[:SCNai].k => 0.0989984,
    sglt[:r12].r => 10.1796,
    sglt[:r23].r => 202.023,
    sglt[:r34].r => 505.058,
    sglt[:r45].r => 8080.93,
    sglt[:r56].r => 67.1184,
    sglt[:r61].r => 8.67804,
    sglt[:r25].r => 0.0061077
    ]
u0 = [
    sglt[:CNao].q => 2.0,
    # sglt[:So].q => 0.5
]
prob = ODEProblem(model, u0, (0.0, 70.0), p, jac = true, sparse = true)
sol = solve(prob)
plot(sol)