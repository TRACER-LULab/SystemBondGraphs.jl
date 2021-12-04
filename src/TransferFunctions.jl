"""
Construct the transfer function for the Bond Graph with `s` as the Laplace Variable for the state-space representation of a linear system of the form 
    ``\\dot{\\vec{x}} = \\boldsymbol{A}\\vec{x}+\\boldsymbol{B}\\vec{u}``   
    ``\\vec{y} = \\boldsymbol{C}\\vec{x} + \\boldsymbol{D}\\vec{u}``   
"""
function state_space(BG::AbstractBondGraph; ps = Dict{Any, Any}())
    if length(states(BG.model)) == 0
        error("Model not generated. Run generate_model!")
    end
    state_vars = states(BG.model)
    eqns = equations(BG.model)
    eqns = map(eqn -> expand(substitute(eqn, ps)), eqns)
    eqns_dict = Dict(get_variables(getfield(eqn, :lhs))[1] => eqn for eqn in eqns)
    obs = observed(BG.model)
    obs = map(eqn -> expand(substitute(eqn, ps)), obs)
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
            A[i,j] = expand(substitute(substitute(eqns_dict[state_vars[i]].rhs, state_dict), in_dict))
            state_dict[state_vars[j]] = 0.0
        end
        for j in eachindex(ins)
            in_dict[ins[j]] = 1.0
            B[i,j] = expand(substitute(substitute(eqns_dict[state_vars[i]].rhs, state_dict), in_dict))
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
            D[i,j] = substitute(substitute(main_expr, state_dict), in_dict) |> simplify
            in_dict[ins[j]] = 0.0
        end
    end
    q = Dict(map(i-> state_vars[i] => i, eachindex(state_vars)))
    w = Dict(map(i-> ins[i] => i, eachindex(ins)))
    e = Dict(map(i-> obs_vars[i] => i, eachindex(obs_vars)))
    return A, B, C, D, q, w, e
end