using Base: OrderStyle
vin(t) = ((U/d)*t<=1) ? (h/d)*π*U*cos(π*U/d*t) : 0.0
ϕ9(q) = (q <= qs0) ? ks1*q : ks1*qs0+ks2*(q-qs0)
ϕ2(q) = (q>=0.0) ? q*kt : 0.0
ϕ8(f) = B*f^3
@register vin(t) 
@register ϕ9(q)
@register ϕ2(q)
@register ϕ8(f)

@variables t e1(t) e2(t) e3(t) e4(t) e5(t) e6(t) e7(t) e8(t) e9(t) e10(t) e11(t) e12(t)
@variables f1(t) f2(t) f3(t) f4(t) f5(t) f6(t) f7(t) f8(t) f9(t) f10(t) f11(t) f12(t)
@variables q2(t) q9(t) p5(t) p12(t)

D = Differential(t)
U = 0.9
d = 1.0
h = 0.25
# Masses
ms = 320
mus = ms/6
# Suspension
fs = 1.0
ωs = 2*π*fs
ks1 = ms*ωs^2
ks2 = 10*ks1
g = 9.81
qs_init = ms*g/ks1
qus_init = (ms+mus)*g/kt
qs0 = 1.3*qs_init
# Tire
kt = 10*ks1
# Damper
B = 1500
eqns = [
    0 ~ f1 - vin(t),
    0 ~ e4 - mus*g,
    0 ~ e11 - ms*g,
    0 ~ e2 - e1,
    0 ~ e2 - e3,
    0 ~ f1 - f3 - f2,
    D(q2) ~ f2,
    0 ~ f5 - f3,
    0 ~ f5 - f4,
    0 ~ f5 - f6,
    0 ~ e3 - e4 - e5 - e6,
    D(p5) ~ e5,
    0 ~ e7 - e6,
    0 ~ e7 - e10,
    0 ~ f6 - f10 - f7,
    0 ~ f7 - f8,
    0 ~ f7 - f9,
    0 ~ e8 + e9 - e7,
    D(q9) ~ f9,
    0 ~ e8 - ϕ8(f8),
    0 ~ f12 - f10,
    0 ~ f12 - f11,
    0 ~ e10 - e11 - e12,
    D(p12) ~ e12,
    0 ~ e2 - ϕ2(q2),
    0 ~ f5 - p5/mus,
    0 ~ e9 - ϕ9(q9),
    0 ~ f12 - p12/ms
]

sys = ODESystem(eqns, t)
sys = structural_simplify(sys)

# Set initial conditions
u0 = [
    q9  => qs_init,
    q2  => qus_init,
    p5 => 0.0,
    p12  => 0.0
]
# Set TimeSpan
tspan = (0.0, 2.0)
simplify_model!(quarterCar)
prob = ODAEProblem(sys, u0, tspan, sparse =true, jac = true)
sol = solve(prob)
Plots.plot(sol, vars = [e2, e9, e8])