"""
Construct the transfer function for the Bond Graph with `s` as the Laplace Variable for a linear system of the form \$\\dot{\\vec{x}} = \\boldsymbol{A}\\vec{x}+\\boldsymbol{B}\\vec{u}\$
"""
function state_space(BG::BondGraph, input, output; ps = Dict{Any, Any}())
    sts = states(BG.model)
    eqns = equations(BG.model)
    eqns = map(eqn -> expand(substitute(eqn, ps)), eqns)
    obs = observed(BG.model)
    obs = map(eqn -> expand(substitute(eqn, ps)), obs)
    obs_dict = Dict(obs[i].lhs => obs[i].rhs for i in eachindex(obs))
    Se_vertices = filter_vertices(BG.graph, :type, :Se)
    Se_sys = map(v -> get_prop(BG.graph, v, :sys).Se, Se_vertices)
    Sf_vertices = filter_vertices(BG.graph, :type, :Sf)
    Sf_sys = map(v -> get_prop(BG.graph, v, :sys).Sf, Sf_vertices)
    ins = reduce(vcat, [Se_sys; Sf_sys])

    in_dict = Dict(ins .=> 0.0)
    st_dict = Dict(sts .=> 0.0)
    A = []
    B = []
    # eqns = map(eqn -> substitute(eqn, ps), eqns)
    for i ∈ eachindex(eqns)
        Aᵢ = []
        for j ∈ eachindex(sts)
            st_dict[sts[j]] = 1.0
            push!(Aᵢ, substitute(substitute(eqns[i].rhs, st_dict), in_dict) |> simplify)
            st_dict[sts[j]] = 0.0
        end
        Bᵢ = []
        for j ∈ eachindex(ins)
            in_dict[ins[j]] = 1.0
            push!(Bᵢ, substitute(substitute(eqns[i].rhs, st_dict), in_dict) |> simplify)
            in_dict[ins[j]] = 0.0
        end
        push!(A, Aᵢ) 
        push!(B, Bᵢ)
    end

    A = reduce(hcat, A) |> permutedims
    B = reduce(hcat, B) |> permutedims
    C = []
    D = []
    main_expr = obs_dict[output]
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
    for j ∈ eachindex(sts)
        st_dict[sts[j]] = 1.0
        push!(C, substitute(substitute(main_expr, st_dict), in_dict) |> expand)
        st_dict[sts[j]] = 0.0
    end
    for j ∈ eachindex(ins)
        in_dict[ins[j]] = 1.0
        push!(D, substitute(substitute(main_expr, st_dict), in_dict) |> expand)
        in_dict[ins[j]] = 0.0
    end
    # TF = (s * I(length(sts)) - A)^(-1) * B
    # println(size(TF))
    # res = Dict()
    # for i ∈ eachindex(sts)
    #     for j ∈ eachindex(ins)
    #         res[sts[i], ins[j]] = TF[i, j]
    #     end
    # end
    # return res
    return A, B, C, D, sts, ins
end

function state_matrices(BG::BondGraph, s; ps = Dict{Any, Any}())
    sts = states(BG.model)
    eqns = equations(BG.model)
    eqns = map(eqn -> simplify(expand(substitute(eqn, ps))), eqns)
    Se_vertices = filter_vertices(BG.graph, :type, :Se)
    Se_sys = map(v -> get_prop(BG.graph, v, :sys).Se, Se_vertices)
    Sf_vertices = filter_vertices(BG.graph, :type, :Sf)
    Sf_sys = map(v -> get_prop(BG.graph, v, :sys).Sf, Sf_vertices)
    ins = reduce(vcat, [Se_sys; Sf_sys])

    in_dict = Dict(ins .=> 0.0)
    st_dict = Dict(sts .=> 0.0)
    A = []
    B = []
    # eqns = map(eqn -> substitute(eqn, ps), eqns)
    for i ∈ eachindex(eqns)
        Aᵢ = []
        for j ∈ eachindex(sts)
            st_dict[sts[j]] = 1.0
            push!(Aᵢ, substitute(substitute(eqns[i].rhs, st_dict), in_dict) |> expand)
            st_dict[sts[j]] = 0.0
        end
        Bᵢ = []
        for j ∈ eachindex(ins)
            in_dict[ins[j]] = 1.0
            push!(Bᵢ, substitute(substitute(eqns[i].rhs, st_dict), in_dict) |> expand)
            in_dict[ins[j]] = 0.0
        end
        push!(A, Aᵢ) 
        push!(B, Bᵢ)
    end

    A = reduce(hcat, A) |> permutedims
    B = reduce(hcat, B) |> permutedims
    sts = Dict([sts[i] => i for i in eachindex(sts)])
    ins = Dict([ins[i] => i for i in eachindex(ins)])
    # A = simplify.(expand.(A))
    # B = simplify.(expand.(B))
    return A, B, sts, ins
end