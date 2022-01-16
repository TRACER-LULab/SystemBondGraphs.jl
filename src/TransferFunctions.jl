"""
Construct the transfer function for the Bond Graph with `s` as the Laplace Variable for the state-space representation of a linear system of the form 
    ``\\dot{\\vec{x}} = \\boldsymbol{A}\\vec{x}+\\boldsymbol{B}\\vec{u}``   
    ``\\vec{y} = \\boldsymbol{C}\\vec{x} + \\boldsymbol{D}\\vec{u}``   
"""
function state_space(BG::AbstractBondGraph, model::ODESystem; ps = Dict{Any,Any}())
    if length(filter(eq -> eq.lhs isa Int64, full_equations(model))) > 0
        de_model, alg_eqns = remove_algebraic(BG, model)
        state_vars = states(de_model)
        eqns = full_equations(de_model)
        eqns = map(eqn -> expand(substitute(eqn, ps)), eqns)
        eqns_dict = Dict(get_variables(getfield(eqn, :lhs))[1] => eqn for eqn in eqns)
        obs = [observed(model); alg_eqns]
        obs = map(eqn -> expand(substitute(eqn, ps)), obs)
    else
        state_vars = states(model)
        eqns = full_equations(model)
        eqns = map(eqn -> expand(substitute(eqn, ps)), eqns)
        eqns_dict = Dict(get_variables(getfield(eqn, :lhs))[1] => eqn for eqn in eqns)
        obs = observed(model)
        obs = map(eqn -> expand(substitute(eqn, ps)), obs)
    end

    obs_dict = Dict(obs[i].lhs => obs[i].rhs for i in eachindex(obs))
    obs_vars = collect(keys(obs_dict))
    Se_vertices = filter_vertices(BG.graph, :type, :Se)
    Se_sys = map(v -> get_prop(BG.graph, v, :sys).Se, Se_vertices)
    Sf_vertices = filter_vertices(BG.graph, :type, :Sf)
    Sf_sys = map(v -> get_prop(BG.graph, v, :sys).Sf, Sf_vertices)
    ins = reduce(vcat, [Se_sys; Sf_sys])
    in_dict = Dict(ins .=> 0.0)
    state_dict = Dict(state_vars .=> 0.0)
    A = Matrix{Num}(undef, length(state_vars), length(state_vars))
    B = Matrix{Num}(undef, length(state_vars), length(ins))
    C = Matrix{Num}(undef, length(obs_vars), length(state_vars))
    D = Matrix{Num}(undef, length(obs_vars), length(ins))
    for i in eachindex(state_vars)
        for j in eachindex(state_vars)
            state_dict[state_vars[j]] = 1.0
            A[i, j] = simplify(substitute(substitute(eqns_dict[state_vars[i]].rhs, state_dict), in_dict)) |> expand
            state_dict[state_vars[j]] = 0.0
        end
        for j in eachindex(ins)
            in_dict[ins[j]] = 1.0
            B[i, j] = expand(simplify(substitute(substitute(eqns_dict[state_vars[i]].rhs, state_dict), in_dict)))
            in_dict[ins[j]] = 0.0
        end
    end

    for i in eachindex(obs_vars)
        main_expr = obs_dict[obs_vars[i]]
        vars = get_variables(main_expr)
        finished = false
        while !finished
            status = 0
            for var ∈ vars
                try
                    main_expr = substitute(main_expr, Dict(var => obs_dict[var]))
                catch
                    status += 1
                end
            end
            finished = (status == length(vars))
            vars = get_variables(main_expr)
        end
        for j ∈ eachindex(state_vars)
            state_dict[state_vars[j]] = 1.0
            C[i, j] = substitute(substitute(main_expr, state_dict), in_dict) |> simplify
            state_dict[state_vars[j]] = 0.0
        end
        for j ∈ eachindex(ins)
            in_dict[ins[j]] = 1.0
            D[i, j] = substitute(substitute(main_expr, state_dict), in_dict) |> simplify
            in_dict[ins[j]] = 0.0
        end
    end
    state_dict = Dict(map(i -> state_vars[i] => i, eachindex(state_vars)))
    in_dict = Dict(map(i -> ins[i] => i, eachindex(ins)))
    obs_dict = Dict(map(i -> obs_vars[i] => i, eachindex(obs_vars)))
    return A, B, C, D, state_dict, in_dict, obs_dict
end