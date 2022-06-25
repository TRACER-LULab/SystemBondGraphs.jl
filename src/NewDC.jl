function resolve_equations(constituitive_equation::Equation, derivative_form::Equation, algebraic_eqns::Vector{Equation}, differential_eqns::Vector{Equation}, ps, svs)
    alg_eqs = algebraic_eqns
    check(x) = (isa(x, Term) && (indexin([x], ps)[1] === nothing) && (indexin([x], svs)[1] === nothing))
    terms = filter(check, get_variables(constituitive_equation.rhs))
    terms_old = []
    while !check(constituitive_equation) && length(intersect(terms, terms_old)) != length(terms) && length(terms) != 0
        terms_old = terms
        for term in terms
            for index in eachindex(alg_eqs)
                eq = 0 ~ alg_eqs[index].lhs - alg_eqs[index].rhs
                eq_terms = filter(x -> isa(x, Term), get_variables(eq))
                found = !(indexin([term], eq_terms)[1] === nothing)
                if found
                    eq = flatten_fractions(expand(eq.rhs))
                    if hasfield(typeof(eq), :num)
                        eq = 0 ~ eq.num
                    else
                        eq = 0 ~ eq
                    end
                    # @show eq, term
                    rhs = Symbolics.solve_for(eq, term)
                    lhs = term
                    constituitive_equation = flatten_fractions(simplify(substitute(constituitive_equation, Dict(lhs => rhs))))
                    # differential_eqns = simplify.(substitute.(differential_eqns, (Dict(lhs => rhs),)))
                    alg_eqs = simplify.(substitute.(alg_eqs, (Dict(lhs => rhs),)))
                    alg_eqs[index] = 0 ~ 0
                    # popat!(alg_eqs, index)
                    # break
                end
            end
        end
        terms = filter(check, get_variables(constituitive_equation.rhs))
    end

    # ff(x) = SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(SymbolicUtils.add_with_div))(x)
    # Fixpoint(Postwalk(add_with_div))(x)
    constituitive_equation = constituitive_equation.lhs ~ flatten_fractions(simplify(constituitive_equation.rhs))
    for index in eachindex(differential_eqns)
        deq = differential_eqns[index]
        sub_dict = Dict(constituitive_equation.lhs => flatten_fractions(simplify(constituitive_equation.rhs)))
        derived_eqn_rhs = flatten_fractions(simplify(substitute(deq.lhs, sub_dict)))
        derived_eqn_rhs = expand_derivatives(derived_eqn_rhs)
        dsub_dict = Dict(map(x -> x.lhs => x.rhs, differential_eqns))
        derived_eqn_rhs = flatten_fractions(simplify(substitute(derived_eqn_rhs, dsub_dict)))
        differential_eqns[index] = deq.lhs ~ derived_eqn_rhs
        dsub = Dict(derivative_form.rhs => derivative_form.lhs)
        differential_eqns[index] = flatten_fractions(simplify(expand(substitute(differential_eqns[index], dsub))))
    end
    # display(differential_eqns)
    # [differential_eqns; constituitive_equation; derived_eqn; algebraic_eqns]
    differential_eqns#; constituitive_equation; algebraic_eqns]

end
function derivative_casuality(BG; skip = [])
    model = BondGraphs.graph_to_model(BG)
    eqs = full_equations(model)
    ps = [parameters(model); skip]
    fn(g, v) = (get_prop(g, v, :type) == :C || get_prop(g, v, :type) == :I) && !get_prop(g, v, :causality)
    svs = map(v -> get_prop(BG.graph, v, :state_var), filter_vertices(BG.graph, fn))
    svs = reduce(vcat, svs)
    dc_nodes = filter_vertices(BG.graph, :causality, :true)
    D = Differential(model.iv)
    differential_eqns, algebraic_eqns = BondGraphs.get_diff_and_alg_eqns(eqs)
    for node in dc_nodes
        name = get_prop(BG.graph, node, :name)
        if get_prop(BG.graph, node, :type) == :C
            constituitive_equation = BG[name].q ~ BG[name].C * BG[name].e
            starting_equation = BG[name].e ~ BG[name].q / BG[name].C
            derivative_form = BG[name].f ~ D(BG[name].q)
        elseif get_prop(BG.graph, node, :type) == :I
            constituitive_equation = BG[name].p ~ BG[name].I * BG[name].f
            starting_equation = BG[name].f ~ BG[name].p / BG[name].I
            derivative_form = BG[name].e ~ D(BG[name].p)
        end
        start_index = indexin([starting_equation], algebraic_eqns)[1]
        filtered_eqs = [algebraic_eqns[1:(start_index-1)]; algebraic_eqns[(start_index+1):end]]
        # display("---------")
        # display(filtered_eqs)
        # display("---------")
        differential_eqns = resolve_equations(constituitive_equation, derivative_form, filtered_eqs, differential_eqns, ps, svs) .|> expand .|> simplify .|> flatten_fractions
        # eqs = [differential_eqns...; algebraic_eqns; [starting_equation]]
        # display(eqs)
    end
    eqs = [differential_eqns; algebraic_eqns]
    model = ODESystem(eqs, model.iv, name = nameof(model))
    return model
end
