### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ ae385a3a-eb33-11eb-1706-b5809ce98d3e
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 50e0db66-0776-4e3e-a907-18fd4f8985f6
begin
	using BondGraphs
	using ModelingToolkit
	using DifferentialEquations
end

# ╔═╡ 499ece3c-d412-4995-a02b-72839e6cb035
md"""
# 1. The Creation of a BondGraph
"""

# ╔═╡ c3bed75a-8cc3-4096-95c6-62743f71b1b3
md"""### Create the independent variable"""

# ╔═╡ b10d3582-b3b8-420a-9d4d-97486a730558
md"""### Create BondGraph"""

# ╔═╡ e0396eee-5ee2-4767-800c-7e499841aa70
@variables t;

# ╔═╡ 05749e70-0fd3-41c7-92fe-ea50d84a8958
visco = BondGraph(t)

# ╔═╡ be3ebd1a-30cc-4763-9100-9e4000a4cd92
begin
	# Inputs
	add_Se!(visco, :σ₁)
	add_Se!(visco, :σ₂)
	add_Se!(visco, :εE)
	add_Bond!(visco, :b2)
	add_Bond!(visco, :b5)
	add_Bond!(visco, :b6)
	add_Bond!(visco, :b8)
	add_Bond!(visco, :b9)
	add_Bond!(visco, :b11)
	## Add Dampers
	add_R!(visco, :R1)
	add_R!(visco, :R2)
end

# ╔═╡ c37a36ca-9777-43d9-8428-bb9caa19784a
md"""### Create the necessary bonds"""

# ╔═╡ 93120f7c-3f6e-4d63-aff8-88bb37b8f2d8
md"""### Create the hyperelastic relationship
"""

# ╔═╡ 415ed7bd-8fd2-48e2-be7b-3f5e4ae1c6d5
function ϕi(𝐞, 𝐪, params)
    λ₁, λ₂ = 𝐪
    λ₃ = 1 / λ₁ / λ₂
    μ = params
	σ1 = μ * (λ₁^2 - λ₃^2)
    σ2 = μ * (λ₂^2 - λ₃^2)
    return [σ1; σ2]
end 

# ╔═╡ 1c33e97e-3d26-4652-8409-02bcb0afede4
md"""### Create Multi-Port Elements"""

# ╔═╡ 36890040-6dd5-4b6e-a7b5-15087fe42869
add_1J!(visco, Dict([
    :σ₁ => false, 
    :εE => false,
    :b5 => true,
    :b2 => true
    ]), :J1_1)

# ╔═╡ 1b73dc3e-1767-441d-ba76-cd9fd1d0e7d1
add_1J!(visco, Dict([
    :σ₂ => false, 
    :εE => false,
    :b8 => true,
    :b11 => true
    ]), :J1_1)

# ╔═╡ 1819df2e-b4a9-4a93-a7db-11d7dd92d4e6
add_0J!(visco, Dict([
    :b2 => false,
    :R1 => true,
    :b6 => true
    ]), :J0_1)

# ╔═╡ eab485fc-eac1-4f83-aaa0-79b4803167bc
add_0J!(visco, Dict([
    :b11 => false,
    :R2 => true,
    :b9 => true
    ]), :J0_1)

# ╔═╡ 08390631-5454-4f21-afe6-b68fc7572892
generate_model!(visco)

# ╔═╡ 08d48fb7-174b-49de-a5a0-109e3d82905c
equations(visco.model)

# ╔═╡ 26ae82cd-8073-4c10-972d-f2cf31865755


# ╔═╡ Cell order:
# ╠═ae385a3a-eb33-11eb-1706-b5809ce98d3e
# ╠═50e0db66-0776-4e3e-a907-18fd4f8985f6
# ╟─499ece3c-d412-4995-a02b-72839e6cb035
# ╠═be3ebd1a-30cc-4763-9100-9e4000a4cd92
# ╟─c3bed75a-8cc3-4096-95c6-62743f71b1b3
# ╠═05749e70-0fd3-41c7-92fe-ea50d84a8958
# ╟─b10d3582-b3b8-420a-9d4d-97486a730558
# ╠═e0396eee-5ee2-4767-800c-7e499841aa70
# ╟─c37a36ca-9777-43d9-8428-bb9caa19784a
# ╟─93120f7c-3f6e-4d63-aff8-88bb37b8f2d8
# ╠═415ed7bd-8fd2-48e2-be7b-3f5e4ae1c6d5
# ╟─1c33e97e-3d26-4652-8409-02bcb0afede4
# ╠═36890040-6dd5-4b6e-a7b5-15087fe42869
# ╠═1b73dc3e-1767-441d-ba76-cd9fd1d0e7d1
# ╠═1819df2e-b4a9-4a93-a7db-11d7dd92d4e6
# ╠═eab485fc-eac1-4f83-aaa0-79b4803167bc
# ╠═08390631-5454-4f21-afe6-b68fc7572892
# ╠═08d48fb7-174b-49de-a5a0-109e3d82905c
# ╠═26ae82cd-8073-4c10-972d-f2cf31865755
