using BondGraphs
using ModelingToolkit
using DifferentialEquations
# using Plots
using LinearAlgebra
## Analysis for Third Damper Model
@parameters t
# Create Empty Bondgraph
third_damper = BondGraph(t)
@parameters α ω
# vin(e, f, t, p) = p[1]*sin(ω*t)
# add_Sf!(third_damper, vin, [α, ω], :vin)
add_Sf!(third_damper, :vin)
add_C!(third_damper, :kt1)
add_Bond!(third_damper, :b3)
add_0J!(third_damper, 
Dict(
    :vin => true,
    :kt1 => false,
    :b3 => false)
    ,
    :J01)
add_I!(third_damper, :mus1)
add_Se!(third_damper, :musg1)
add_Bond!(third_damper, :b6)
add_Bond!(third_damper, :b24)
add_Bond!(third_damper, :b23)
add_1J!(third_damper, Dict(
    :b3 => true, 
    :mus1 => false,
    :musg1 => true,
    :b6 => false,
    :b23 => false,
    :b24 => false
    ),
    :J11)
add_Bond!(third_damper, :b7)
add_Bond!(third_damper, :b10)
add_0J!(third_damper, Dict(
    :b6 => true, 
    :b7 => false,
    :b10 => false
    ),
    :J02)
add_C!(third_damper, :k1)
add_R!(third_damper, :b1)
add_1J!(third_damper, Dict(
    :k1 => false,
    :b1 => false,
    :b7 => true
    ), 
    :J12)
add_Se!(third_damper, :msg)
add_I!(third_damper, :ms)
add_Bond!(third_damper, :b13)
add_1J!(third_damper, Dict(
    :b10 => true, 
    :msg => false,
    :ms => false,
    :b13 => true
    ),
    :J13)
add_Bond!(third_damper, :b14)
add_Bond!(third_damper, :b17)
add_0J!(third_damper, Dict(
    :b13 => false,
    :b14 => false,
    :b17 => true
    ),
    :J03)
add_C!(third_damper, :k2)
add_R!(third_damper, :b2)
add_1J!(third_damper, Dict(
    :b14 => true, 
    :k2 => false,
    :b2 => false
    ), 
    :J14)
add_Se!(third_damper, :musg2)
add_I!(third_damper, :mus2)
add_C!(third_damper, :kt2)
add_Bond!(third_damper, :b21)
add_Bond!(third_damper, :b27)
add_1J!(third_damper, Dict(
    :b17 => false,
    :musg2 => true,
    :mus2 => false,
    :kt2 => false,
    :b21 => true,
    :b27 => true
    ), 
    :J15)
add_R!(third_damper, :b)
add_0J!(third_damper, Dict(
    :b24 => true, 
    :b => false,
    :b27 => false
    ), 
    :J04)
add_C!(third_damper, :k)
add_0J!(third_damper, Dict(
    :b23 => true, 
    :b21 => false,
    :k => false
    ),
    :J05)

## Generate and Simplify Model
generate_model!(third_damper)
third_damper.model = structural_simplify(third_damper.model)

# Set Parameters for study
@variables mus Sf1
ps = Dict{Num , Real}(
    # third_damper[:vin].α   => 0.01,
    # third_damper[:vin].ω   => 1.0,
    third_damper[:ms].I     => 680.0,
    third_damper[:mus1].I   => mus,
    third_damper[:mus2].I   => mus,
    third_damper[:musg1].Se   => mus*9.81,
    third_damper[:musg2].Se   => mus*9.81,
    third_damper[:msg].Se   => 680.0*9.81,
    third_damper[:k1].C     => 1/32_000,
    third_damper[:k2].C     => 1/32_000,
    third_damper[:kt1].C    => 1/360_000,
    third_damper[:kt2].C    => 1/360_000,
    third_damper[:k].C      => 1/360_000,
    third_damper[:b1].R     => 2798.86,
    third_damper[:b2].R     => 2798.86,
    third_damper[:b].R      => 1119.54,
    )

u0 = Dict{Num, Real}(
    third_damper[:k].q => 0.0,
    third_damper[:k1].q => 0.0,
    third_damper[:k2].q => 0.0,
    third_damper[:kt1].q => 0.0,
    third_damper[:kt2].q => 0.0,
    third_damper[:mus1].p => 0.0,
    third_damper[:mus2].p => 0.0,
    third_damper[:ms].p => 0.0,#
)

tspan = (0.0, 20.0)
prob = ODAEProblem(third_damper.model, u0, tspan, ps)
sol = solve(prob)

## Get A & B state_matrices
A, B, C, D, sts, ins, obs = state_space(third_damper, ps = ps)
_A_func = eval(build_function(A,[mus],parallel=Symbolics.MultithreadedForm())[1])
AF(mus) = reshape(_A_func([mus]), size(A))
_B_func = eval(build_function(B,[mus],parallel=Symbolics.MultithreadedForm())[1])
BF(mus) = reshape(_B_func([mus]), size(B))
_C_func = eval(build_function(C,[mus],parallel=Symbolics.MultithreadedForm())[1])
CF(mus) = reshape(_C_func([mus]), size(C))
_D_func = eval(build_function(D,[mus],parallel=Symbolics.MultithreadedForm())[1])
DF(mus) = reshape(_D_func([mus]), size(D))
## Adjustment to Match MATLAB Code
# A[2, 6] = -A[2, 6]
# A[6, 2] = -A[6, 2]

## Plotting
AR_plot = plot()
PA_plot = plot()
linetype = [:solid, :dash, :dot]
m = [25.0, 50.0, 75.0]
in_var = third_damper[:vin].Sf
out_var = third_damper[:ms].e
for i in eachindex(m)
    freqs = 10 .^(0:0.01:3)
    _A = AF(m[i])
    _B = BF(m[i])
    _C = CF(m[i])
    _D = DF(m[i])
    TF(s) = (_C*(s*I(size(_A, 1))-_A)^(-1)*_B+_D)[obs[out_var], ins[in_var]]
    res = TF.(freqs.*1im)./ps[third_damper[:b1].R]
    AR = abs.(res)
    PA = rad2deg.(angle.(res))
    plot!(AR_plot, freqs./2π, AR, label = "3rd Damper - "*string(m[i]), linestyle = linetype[i], color = :blue)
    plot!(PA_plot, freqs./2π, PA, label = "3rd Damper - "*string(m[i]), linestyle = linetype[i], color = :blue)
end
# scatter(fft(AR))
# plot!(phase)
plot!(AR_plot, size = (400,100).*2, xtick = 10 .^ (0:1:5), legend = :topleft, xscale = :log10, xlabel = "Frequency [rad/s]", ylabel = "AR [(Fₘₛ+mₛg)/(Vᵢₙb₁)]")
plot!(PA_plot, size = (400,100).*2, xtick = 10 .^ (0:1:5), legend = :bottomleft, xscale = :log10, xlabel = "Frequency [rad/s]", ylabel = "Phase Angle [∘]")
plot!(AR_plot, PA_plot, layout = (2, 1), size = (800, 500))

## Normal Car

@parameters t
car = BondGraph(t)
# Build Bond Graph
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
add_Bond!(car, :b6)
add_Bond!(car, :b23)
add_1J!(car, Dict(
    :b3 => true,
    :mus1 => false,
    :mus1g => false,
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
add_1J!(car, Dict(
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
add_1J!(car, Dict(
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
    :mus2g => false, 
    :mus2 =>false,
    :kt2 => false,
    :b21 => true
    ),
    :J15)
add_C!(car, :k)
add_0J!(car, Dict(
    :b23 => true,
    :k => false,
    :b21 => false
    ),
    :J04) 

## Generate Model
generate_model!(car)
car.model = structural_simplify(car.model)

# Create Parameters
@variables mus
ps = Dict{Num , Real}(
    car[:ms].I     => 680.0,
    car[:mus1].I   => mus,
    car[:mus2].I   => mus,
    car[:k1].C     => 1/32_000,
    car[:k2].C     => 1/32_000,
    car[:kt1].C    => 1/360_000,
    car[:kt2].C    => 1/360_000,
    car[:k].C      => 1/360_000,
    car[:b1].R     => 2798.86,
    car[:b2].R     => 2798.86,
    )

# Get A & B state_matrices
@variables s
A, B, C, D, sts, ins, obs = state_space(car, ps = ps)
_A_func = eval(build_function(A,[mus],parallel=Symbolics.MultithreadedForm())[1])
AF(mus) = reshape(_A_func([mus]), size(A))
_B_func = eval(build_function(B,[mus],parallel=Symbolics.MultithreadedForm())[1])
BF(mus) = reshape(_B_func([mus]), size(B))
_C_func = eval(build_function(C,[mus],parallel=Symbolics.MultithreadedForm())[1])
CF(mus) = reshape(_C_func([mus]), size(C))
_D_func = eval(build_function(D,[mus],parallel=Symbolics.MultithreadedForm())[1])
DF(mus) = reshape(_D_func([mus]), size(D))
# Plot
m = [25.0, 50.0, 75.0]
# AR_plot = plot()
# PA_plot = plot()
##
linetype = [:solid, :dash, :dot]
in_var = car[:vin].Sf
out_var = car[:ms].e
for i ∈ eachindex(m)
    freqs = 10 .^(0:0.01:3)
    _A = AF(m[i])
    _B = BF(m[i])
    _C = CF(m[i])
    _D = DF(m[i])
    TF(s) = (_C*(s*I(size(_A, 1))-_A)^(-1)*_B+_D)[obs[out_var], ins[in_var]]
    res = TF.(freqs.*1im)./ps[car[:b1].R]
    AR = abs.(res)
    PA = rad2deg.(atan.(imag.(res), real.(res)))
    plot!(AR_plot, freqs./2π, AR, label = "Conventional - "*string(m[i]), linestyle = linetype[i], color = :red)
    plot!(PA_plot, freqs./2π, PA, label = "Conventional - "*string(m[i]), linestyle = linetype[i], color = :red)
end

plot!(AR_plot, size = (400,100).*2, xtick = 10 .^ (0:1:5), xscale = :log10)
plot!(PA_plot, size = (400,100).*2, xtick = 10 .^ (0:1:5), xscale = :log10)
plot!(AR_plot, PA_plot, layout = (2, 1), size = (800, 500))