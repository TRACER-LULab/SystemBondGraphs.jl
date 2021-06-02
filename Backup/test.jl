using BondGraph
#########
## Mass Spring Damper systems
@variables t e(t) f(t) q(t) p(t)
msd = BondGraph(t)
# Elements
function rf(f,t)
    Int(floor(t))
end

@register rf(f, t)
add_R!(msd, rf, :r1)
add_C!(msd, 1.0, :c1)
add_I!(msd, 1.0, :i1)
add_Se!(msd, 1.0, :se)
# Junctions
add_1J!(msd, Dict([
    :r1=>true, 
    :c1=>true, 
    :i1=>true, 
    :se=>false
    ]), :J1)
# Creating Model
generate_model!(msd)
simplify_model!(msd)
ivs = get_states(msd)
ic = Dict(map(x->x=>0.0, ivs))
set_conditions!(msd, ic)
prob = generate_ODE(msd, (0.0, 10.0))
sol = solve(prob)
plot(sol)

## ## ## ## ## ##
## Quarter Car ##
## ## ## ## ## ##
# Create Bond Graph
quarterCar = BondGraph(t)
# Add Elements
@parameters ms mus fs ωs ks1 ks2 kt qs_init qus_init qs0 B h d U g
v(t) = ((U/d)*t<=1) ? (h/d)*π*U*cos(π*U/d*t) : 0.0
@register v(t)
add_Sf!(quarterCar, v(t), :Vin)
add_C!(quarterCar, kt, :kt)
add_Bond!(quarterCar, :b3)
add_Se!(quarterCar, sin(t), :Fus)
add_I!(quarterCar, mus, :mus)
add_Bond!(quarterCar, :b6)
add_Bond!(quarterCar, :b7)
add_R!(quarterCar, B, :b)
function ϕ(q, e, t)
    if q < qs0
        return e/ks1
    else
        return (Fs-ks1*qs0)/ks2+qs0
    end
end
@register ϕ(q, e,t)
add_C!(quarterCar, ϕ, :ks)
add_Bond!(quarterCar, :b10)
add_Se!(quarterCar, ms*g, :Fs)
add_I!(quarterCar, ms, :ms)
# Add Junctions
add_0J!(quarterCar, Dict([
    :Vin=>false, 
    :kt=>:true, 
    :b3=>true]),
    :J01)
add_1J!(quarterCar, Dict([
    :b3=>false, 
    :mus=>true, 
    :Fus=>true,
    :b6=>true]), 
    :J11)
add_0J!(quarterCar, Dict([
    :b6=>false, 
    :b7=>true, 
    :b10=>true]), 
    :J02)
add_1J!(quarterCar, Dict([
    :b7=>false, 
    :ks=>true,
    :b=>true]), 
    :J12)
add_1J!(quarterCar, Dict([
    :b10=>false,
    :Fs=>true,
    :ms=>true
    ]),
    :J13)
# Create Model

generate_model!(quarterCar)
simplify_model!(quarterCar)
ivs = get_states(quarterCar)
ics = Dict([
    kt.q=>ms*g/ks1,
    kus.q=>(ms+mus)*g/kt,
    ms.p=>0.0,
    mus.p=>0.0
    
])
ic = Dict(map(x->x=>0.0, ivs))
set_conditions!(quarterCar, ic)
prob = generate_ODE(quarterCar, (0.0, 10.0))
sol = solve(prob, Rodas4())
plot(sol)

## 
fig6_57 = BondGraph(t)
Rw = 1.0
T = 0.5
kτ = 1/10.0
J = 0.0117
# Create Elements
add_Se!(fig6_57, sin(T^2/Rw/J*t), :ec)
add_R!(fig6_57, Rw, :Rw)
add_C!(fig6_57, kτ, :kτ)
add_I!(fig6_57, J, :J)
add_Se!(fig6_57, sin(10t), :τd)
add_Bond!(fig6_57, :b3)
add_Bond!(fig6_57, :b4)
add_Bond!(fig6_57, :b6)
# Create Junctions
add_1J!(fig6_57, Dict([
    :ec => false, 
    :Rw => true,
    :b3 => true
    ]), :J11)
add_GY!(fig6_57, T, Dict([
    :b3 => false,
    :b4 => true
    ]), :GY)
add_0J!(fig6_57, Dict([
    :b4 => false, 
    :kτ => true, 
    :b6 => true
    ]), :J01)
add_1J!(fig6_57, Dict([
    :b6 => false,
    :J => true, 
    :τd => false
    ]), :J12)

generate_model!(fig6_57)
simplify_model!(fig6_57)
ivs = get_states(fig6_57)
ic = Dict(map(x->x=>0.0, ivs))
set_conditions!(fig6_57, ic)
prob = generate_ODE(fig6_57, (0.0, 10.0))
sol = solve(prob, Rodas4())
