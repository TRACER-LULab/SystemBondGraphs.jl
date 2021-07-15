using DrWatson
@quickactivate "Bond Graph Modeling"
DrWatson.greet()
using Pkg
Pkg.activate(".")
Pkg.instantiate()
## Package Imports
using Plots 
using ReinforcementLearning
using StableRNGs
using Random
using Flux
using Flux.Losses
using JLD2
using ModelingToolkit
using Symbolics
using DifferentialEquations
using UnPack
## Load BondGraph Model and create symbolic variables
@unpack model, independent_vars, state_vars, params = load(datadir("sims", "ODEModels", "cart_pole_model.jld2"))
## Recreate Independent Variables 
command = "@variables "*string(map(x->string(x)*" ", independent_vars)...)
eval(Meta.parse(command))
## Recreate State Variables
command = "@variables "*string(map(x->string(x)*" ", state_vars)...)
eval(Meta.parse(command))
## Recreate parameters
command = "@parameters "*string(map(x->string(x)*" ", params)...)
eval(Meta.parse(command))
## 
function createEnv(θ0, mass_pole, mass_cart, pole_length, gravity, dt)
    u0 = [
        J₊p => 0.0,    
        mc₊p => 0.0,
        x => 0.0,
        θ => θ0
        ] |> Dict

    ps = [
        mpx₊I => mass_pole,
        mpy₊I => mass_pole,
        mc₊I => mass_cart,
        J₊I => mass_pole * pole_length^2,
        l => pole_length,
        mpg₊Se => mass_pole * gravity,
        r1₊R => 0.01,
        in₊Se => 0.0
        ] |> Dict
    
    prob = ODAEProblem(model, u0, (0.0, 10.0), ps, dt=dt)
    return init(prob, Tsit5(), reltol=1e-6)
end
## Create RL Components
struct CartPoleEnvParams{T}
    gravity::T
    mass_cart::T
    mass_pole::T
    pole_length::T
    force_mag::T
    dt::T
    theta_threshold::T
    x_threshold::T
    theta_start::T
    max_steps::Int
end

Base.show(io::IO, params::CartPoleEnvParams) = print(
    io,
    join(["$p=$(getfield(params, p))" for p in fieldnames(CartPoleEnvParams)], ","),
)

mutable struct CartPoleEnv{T,R <: AbstractRNG} <: AbstractEnv
    params::CartPoleEnvParams{T}
    state::Array{T,1}
    action::Int
    done::Bool
    t::Int
    rng::R
    de_env::OrdinaryDiffEq.ODEIntegrator
end

function CartPoleEnv(;
    T=Float64,
    gravity=9.8,
    mass_cart=1.0,
    mass_pole=0.1,
    pole_length=0.5,
    force_mag=10.0,
    dt=0.02,
    theta_threshold=3.14159 / 4,
    x_threshold=1.0,
    theta_start=randn() * 5 * π / 180,
    max_steps=200,
    rng=Random.GLOBAL_RNG,
)
    params = CartPoleEnvParams{T}(
        gravity,
        mass_cart,
        mass_pole,
        pole_length,
        force_mag,
        dt,
        theta_threshold,
        x_threshold,
        theta_start,
        max_steps,
    )
    high = cp = CartPoleEnv(params, zeros(T, 4), 2, false, 0, rng, createEnv(theta_start, mass_pole, mass_cart, pole_length, gravity, dt))
    reset!(cp)
    cp
end

CartPoleEnv{T}(; kwargs...) where {T} = CartPoleEnv(; T=T, kwargs...)

function RLBase.reset!(env::CartPoleEnv{T}) where {T <: Number}
    env.state[:] = T.([0.0, 0.0, 0.0, randn() * 5 * π / 180])
    env.t = 0
    env.action = 0
    env.done = false
    reinit!(env.de_env, env.state)
    nothing
end

# RLBase.action_space(env::CartPoleEnv) = Base.OneTo(3)ReinforcementLearningBaseReinforcementLearningBase

RLBase.action_space(env::CartPoleEnv{T}) where {T} = Base.OneTo(3)

RLBase.state_space(env::CartPoleEnv{T}) where {T} = Space(
    ClosedInterval{T}[
        (-env.params.x_threshold)..(env.params.x_threshold),
        -1e38..1e38,
        (-env.params.theta_threshold)..(env.params.theta_threshold),
        -1e38..1e38,
    ],
)

RLBase.reward(env::CartPoleEnv{T}) where {T} = env.done ? zero(T) : one(T)
# RLBase.reward(env::CartPoleEnv{T}) where {T} = (3.141549 - abs(mod(env.de_env.sol[θ][end], 3.141549)))^2
RLBase.is_terminated(env::CartPoleEnv) = env.done
RLBase.state(env::CartPoleEnv) = env.state

function (env::CartPoleEnv)(a)
    @assert 0 < a < 4
    env.action = a
    env.t += 1
    env.de_env.p[4] = range(-1, 1, length=3)[a] * env.params.force_mag
    step!(env.de_env, env.params.dt, true)
    env.state[1] = env.de_env.sol[x][end]
    env.state[2] = env.de_env.sol[mc₊f][end]
    env.state[3] = env.de_env.sol[θ][end] # Theta
    env.state[4] = env.de_env.sol[J₊f][end]
    env.done =
        abs(env.state[1]) > env.params.x_threshold ||
        abs(env.state[3]) > env.params.theta_threshold ||
        env.t > env.params.max_steps
    nothing
end

Random.seed!(env::CartPoleEnv, seed) = Random.seed!(env.rng, seed)

## Train the RL Agent
function RL.Experiment(
    ::Val{:JuliaRL},
    ::Val{:BasicDQN},
    ::Val{:CartPoleBG},
    ::Nothing;
    seed=123,
    )
    rng = StableRNG(seed)
    env = CartPoleEnv(; T=Float64, max_steps=400, dt=0.02, rng=rng)
    ns, na = length(state(env)), length(action_space(env))

    agent = Agent(
        policy=QBasedPolicy(
            learner=BasicDQNLearner(
                approximator=NeuralNetworkApproximator(
                    model=Chain(
                        Dense(ns, 128, relu; init=glorot_uniform(rng)),
                        Dense(128, 128, relu; init=glorot_uniform(rng)),
                        Dense(128, na; init=glorot_uniform(rng))
                    ) |> cpu,
                    optimizer=ADAM()
                ),
                batch_size=32,
                min_replay_history=100,
                loss_func=huber_loss,
                rng=rng
            ),
            explorer=EpsilonGreedyExplorer(
                kind=:exp,
                ϵ_stable=0.01,
                decay_steps=500,
                rng=rng
            ),
        ),
        trajectory=CircularArraySARTTrajectory(
            capacity=1000,
            state=Vector{Float32} => (ns,),
        ),
    )
    stop_condition = StopAfterStep(200, is_show_progress=!haskey(ENV, "CI"))
    hook = TotalRewardPerEpisode()
    Experiment(agent, env, stop_condition, hook, "# BasicDQN > CartPole")
end

ex = E`JuliaRL_BasicDQN_CartPoleBG`
run(ex)
## Save
@tagsave(
    datadir("sims", "NNModel", "RL_BasicDQN_CartPole.bson"),
    Dict("model" => ex.policy.policy)
)
##
plot(ex.hook.rewards)
##


anim = @animate for i ∈ eachindex(ex.env.de_env.sol.t)
    plot(
        [ex.env.de_env.sol[x][i], ex.env.de_env.sol[x][i] - sin(ex.env.de_env.sol[θ][i]) * 0.5],
        [0.0, cos(ex.env.de_env.sol[θ][i]) * 0.5], 
        xlims=(-1.5, 1.5),
        ylims=(-1.5, 1.5),
        title="t=" * string(round(ex.env.de_env.sol.t[i], digits=3)),
        linewidth=3,
        aspect_ratio=1
    )
    scatter!(
        [ex.env.de_env.sol[x][i]],
        [0],
        markersize=6,
        c=:red,
        legend=false
    )
end
gif(anim, "cart_pole_BG_damper.gif", fps=50)