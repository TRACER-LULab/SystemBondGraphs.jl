```
transfer_function(BG::BondGraph, s)
Construct the transfer function for the Bond Graph with `s` as the Laplace Variable
```
function transfer_function(BG::BondGraph, s)
    sts = states(BG.model)
    eqns = equations(BG.model)
    Se_vertices = filter_vertices(BG.graph, :type, :Se)
    Se_sys = map(v -> get_prop(BG.graph, v, :sys), Se_vertices)
    Sf_vertices = filter_vertices(BG.graph, :type, :Sf)
    Sf_sys = map(v -> get_prop(BG.graph, v, :sys), Sf_vertices)
    ins = reduce(vcat, map(x -> parameters(x), [Se_sys; Sf_sys]))

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
    (s * I(length(sts)) - A)^(-1) * B
end