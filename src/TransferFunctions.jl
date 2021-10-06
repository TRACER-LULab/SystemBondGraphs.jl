"""
Construct the transfer function for the Bond Graph with `s` as the Laplace Variable for a linear system of the form \$\\dot{\\vec{x}} = \\boldsymbol{A}\\vec{x}+\\boldsymbol{B}\\vec{u}\$
"""
function transfer_function(BG::BondGraph, s; ps = Dict{Any, Any}())
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
    TF = (s * I(length(sts)) - A)^(-1) * B
    println(size(TF))
    res = Dict()
    for i ∈ eachindex(sts)
        for j ∈ eachindex(ins)
            res[sts[i], ins[j]] = TF[i, j]
        end
    end
    return res
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
    return A, B, sts, ins
end

function TF_AB(A, B, s_val, in_idx, out_idx; ps = Dict())
    _A = similar(A)
    _B = similar(B)
    for i ∈ 1:size(A, 1)
        for j ∈ 1:size(A, 2)
            _A[i, j] = substitute(A[i,j], ps)
        end
        for j ∈ 1:size(B, 2)
            _B[i,j] = substitute(B[i,j], ps)
        end
    end
    res = (s_val*I(size(_A, 1))-_A)^(-1)*(_B)
    res[in_idx, out_idx]
end