using BondGraphs
using ModelingToolkit
using DifferentialEquations
using Plots
using Symbolics
using LinearAlgebra
##
@variables t θ(t) ϕ(t) ψ(t) x(t) y(t) z(t) v_b_x(t) v_b_y(t) v_b_z(t) p(t) q(t) r(t) I_xx I_yy I_zz ω[1:4](t) k l b
D = Differential(t)
𝛇 = [x, y, z]
𝛈 = [ϕ, θ, ψ]
𝐪 = [𝛇, 𝛈]
𝐕_b = [v_b_x, v_b_y, v_b_z]
𝛎 = [p, q, r]
C(x) = cos(x)
S(x) = sin(x)
T(x) = tan(x)
𝐑 = [C(ψ)C(θ) C(ψ)S(θ)S(ϕ) - S(ψ)C(ϕ) C(ψ)S(θ)C(ϕ) + S(ψ)S(ϕ);
     S(ψ)C(θ) S(ψ)S(θ)S(ϕ) + C(ψ)C(ϕ) S(ψ)S(θ)C(ϕ) - C(ψ)S(ϕ);
     -S(θ)    C(θ)S(ϕ)              C(θ)C(ϕ)
    ]
𝐖η = [1 0     -S(θ);
      0 C(ϕ)  C(θ)S(ϕ);
      0 -S(ϕ) C(θ)C(ϕ)
     ]

𝐈 = Diagonal([I_xx, I_yy, I_zz])
T_b = sum(x -> k * x^2, ω)
𝐓 = [0, 0, T_b]
τ_ϕ = l * k * (-ω[2]^2 + ω[4]^2)
τ_θ = l * k * (-ω[1]^2 + ω[3]^2)
τ_ψ = sum(x -> b * x, ω)
𝛕_b = [τ_ϕ, τ_θ, τ_ψ]
𝐆 = [0 0 -9.81]
