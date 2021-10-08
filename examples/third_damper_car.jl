using BondGraphs
using ModelingToolkit
using DifferentialEquations
using Plots

## Problem Independent Variable
@parameters t

## Create Empty Bondgraph
koenig = BondGraph(t)
@parameters ω
add_Sf!(koenig, :vin)
add_C!(koenig, :kt1)
add_Bond!(koenig, :b3)
add_0J!(koenig, Dict(
    :vin => true,
    :kt1 => false,
    :b3 => false),
    :J01)
add_I!(koenig, :mus1)
add_Se!(koenig, :musg1)
add_Bond!(koenig, :b6)
add_Bond!(koenig, :b24)
add_Bond!(koenig, :b23)
add_1J!(koenig, Dict(
    :b3 => true, 
    :mus1 => false,
    :musg1 => true,
    :b6 => false,
    :b23 => false,
    :b24 => false
    ),
    :J11)
add_Bond!(koenig, :b7)
add_Bond!(koenig, :b10)
add_0J!(koenig, Dict(
    :b6 => true, 
    :b7 => false,
    :b10 => false
    ),
    :J02)
add_C!(koenig, :k1)
add_R!(koenig, :b1)
add_1J!(koenig, Dict(
    :k1 => false,
    :b1 => false,
    :b7 => true
    ), 
    :J12)
add_Se!(koenig, :msg)
add_I!(koenig, :ms)
add_Bond!(koenig, :b13)
add_1J!(koenig, Dict(
    :b10 => true, 
    :msg => false,
    :ms => false,
    :b13 => true
    ),
    :J13)
add_Bond!(koenig, :b14)
add_Bond!(koenig, :b17)
add_0J!(koenig, Dict(
    :b13 => false,
    :b14 => false,
    :b17 => true
    ),
    :J03)
add_C!(koenig, :k2)
add_R!(koenig, :b2)
add_1J!(koenig, Dict(
    :b14 => true, 
    :k2 => false,
    :b2 => false
    ), 
    :J14)
add_Se!(koenig, :musg2)
add_I!(koenig, :mus2)
add_C!(koenig, :kt2)
add_Bond!(koenig, :b21)
add_Bond!(koenig, :b27)
add_1J!(koenig, Dict(
    :b17 => false,
    :musg2 => true,
    :mus2 => false,
    :kt2 => false,
    :b21 => true,
    :b27 => true
    ), 
    :J15)
add_R!(koenig, :b)
add_0J!(koenig, Dict(
    :b24 => true, 
    :b => false,
    :b27 => false
    ), 
    :J04)
add_C!(koenig, :k)
add_0J!(koenig, Dict(
    :b23 => true, 
    :b21 => false,
    :k => false
    ),
    :J05)
## Generate and Simplify Model
generate_model!(koenig)
koenig.model = structural_simplify(koenig.model)
## Set Parameters for study
@variables mus
ps = Dict{Num , Real}(
    koenig[:ms].I     => 680.0,
    koenig[:msg].Se   => 680.0*9.81,
    koenig[:mus1].I   => mus,
    koenig[:mus2].I   => mus,
    koenig[:musg1].Se => mus*9.81,
    koenig[:musg2].Se => mus*9.81,
    koenig[:k1].C     => 1/32_000,
    koenig[:k2].C     => 1/32_000,
    koenig[:kt1].C    => 1/360_000,
    koenig[:kt2].C    => 1/360_000,
    koenig[:k].C      => 1/360_000,
    koenig[:b1].R     => 2798.86,
    koenig[:b2].R     => 2798.86,
    koenig[:b].R      => 1119.54,
    )

## Create Transfer Function
@variables s
tf1 = transfer_function(koenig, s, ps = ps)[koenig[:ms].p, koenig[:vin].Sf];
tf2 = transfer_function(koenig, s, ps = ps)[koenig[:ms].p, koenig[:msg].Se];
get_variables(tf1)
get_variables(tf2)

## TF = (P_ms*s/V_in + (Se_ms/P_ms)*(P_ms)/(V_in))*(1/b₁)
TF = (tf1*s+tf1/tf2)/ps[koenig[:b1].R];
# TF = (tf1*s);
get_variables(TF)
## Use substitute instead of build_function to get results
function bypass(tf, mus_val, s_val)
    @variables mus, s
    d = Dict(
        mus => mus_val,
        s => s_val
    )
    substitute(tf, d)
end
ans = abs(bypass(TF, 25.0, 10*2*π*1im))
##
AR_plt = plot()
PA_plt = plot()
freqs = 10 .^(0:0.05:3)
for m ∈ [25.0, 50.0, 75.0]
    res = map(f-> bypass(TF, m, f*(2π)*1.0im), freqs)
    display(m)
    AR = getfield.(abs.(res), :val)
    PA = getfield.(angle.(res).*180/π, :val)
    plot!(AR_plt, freqs, AR,  label = string(m))
    display(AR_plt)
    plot!(PA_plt, freqs, PA, label = string(m))
    display(PA_plt)
end
plot!(AR_plt, size = (400,100).*2, xtick = 10 .^ (0:1:3), xscale = :log10, ylims =(0, 4))
plot!(PA_plt, size = (400,100).*2, xtick = 10 .^ (0:1:3), xscale = :log10, ylims = (-400,100))
plot!(AR_plt, PA_plt, layout = (2, 1), size = (800, 400))
##
car = BondGraph(t)

add_Sf!(car, :vin)
add_C!(car, :kt1)
add_Bond!(car, :b3)
add_0J!(car, Dict(
    :vin => true,
    :kt1 => false,
    :b3  => false
    ),
    :J01)
add_I!(car, :mus1)
add_Se!(car, :mus1g)
add_Bond(car, :b6)
add_Bond!(car, :b23)
add_1J!(car, Dict(
    :b3 => true,
    :mus1 => false,
    :mus1g => true,
    :b6 => false,
    :b23 => false
    ),
    :J11)
add_Bond!(car, :b7)
add_Bond!(car, :b10)
add_0J!(car, Dict(
    :b6 => true, 
    :b7 => false, 
    :b10 => false
    ), 
    :J02)
add_C!(car, :k1)
add_R!(car, :b1)
add_1J!(car, DicT(
    :b7 => true, 
    :k1 => false,
    :b1 => false
    ),
    :J12)
add_Se!(car, :msg)
add_I!(car, :ms)
add_Bond!(car, :b13)
add_1J!(car, Dict(
    :b10 => true, 
    :b13 => true,
    :msg => false,
    :ms  => false
    ), 
    :J13)
add_Bond!(car, :b14)
add_Bond!(car, :b17)
add_0J!(car, Dict(
    :b13 => false,
    :b14 => false,
    :b17 => true
    ), 
    :J03)
add_C!(car, :k2)
add_R!(car, :b2)
add_1j!(car, Dict(
    :b14 => true,
    :k2 => false,
    :b2 => false
    ), 
    :J14)
add_Se!(car, :mus2g)
add_I!(car, :mus2)
add_C!(car, :kt2)
add_Bond!(car, :b21)
add_1J!(car, Dict(
    :b17 => false,
    :mus2g => true, 
    :mus2 =>false,
    :kt2 => false,
    :b21 => true
    ),
    :J15)
add_C!(car, :k)
add_0J!(car, DicT(
    :b23 => true,
    :k => false,
    :b21 => false
    ),
    :J04)

generate_model!(car)
car.model = structural_simplify(car)