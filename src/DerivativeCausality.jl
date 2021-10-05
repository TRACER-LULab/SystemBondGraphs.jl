## Resolve Implicit Equations and Derivative causality
function get_args(term)
    # display("start")
    arg_sub = SymbolicUtils.arguments(term)
    args = []
    for arg ∈ arg_sub
        if !(arg isa Number) && !(arg isa Sym)
            if length(SymbolicUtils.arguments(arg)) > 1
                push!(args, get_args(arg)...)
            else 
                push!(args, arg)
            end
        end
    end
    # display("end")
    return args
end

function check_lhs(eqn, term)
    # LHS Check
    if eqn.lhs isa Term
        if (eqn.lhs - term == 0) isa Bool
            return eqn.lhs ~ eqn.rhs
        elseif (eqn.lhs + term == 0) isa Bool
            return eqn.lhs ~ -1 * eqn.rhs
        end
    end
    return false
end

function check_rhs(eqn, term)
    for arg ∈ get_variables(eqn.rhs)
        if (((arg - term == 0) isa Bool) || ((arg + term == 0) isa Bool))
            return term ~ Symbolics.solve_for(eqn, term)
        end
    end
    return false 
end

function find_eqn(term, eqns)
    eqn = Equation[]
    for i ∈ eachindex(eqns)
        res_lhs = check_lhs(eqns[i], term)
        if res_lhs isa Bool
            res_rhs = check_rhs(eqns[i], term)
            if !(res_rhs isa Bool)
                sub_dict[res_rhs.lhs] = res_rhs.rhs
                break
            end
        else 
            return res_lhs.lhs => res_lhs.rhs
        end
    end
end

function get_diff(eqns)
    diff_eqns = []
    alg_eqns = []
    for eqn ∈ eqns
        if !(eqn.lhs isa Real)
            # Check for differential equations
            if SymbolicUtils.operation(eqn.lhs)  isa Differential
                push!(diff_eqns, eqn)
            # Check for equations with state variables
            else 
                push!(alg_eqns, 0 ~ eqn.lhs - eqn.rhs)
            end
        else 
            push!(alg_eqns, eqn)
        end
    end
    return diff_eqns, alg_eqns
end

function get_implicit(state_vars, srch_eqns, alg_eqns)
    res_eqns = []
    for eqn ∈ srch_eqns
        sub_dict = Dict([])
        old_len = -1.0
        eqns = copy(alg_eqns)
        search_terms = [get_args(eqn.rhs)[1]]
        while length(keys(sub_dict)) != old_len
            old_len = length(keys(sub_dict))
            for search_term ∈ search_terms
                for i ∈ eachindex(eqns)
                    res_lhs = check_lhs(eqns[i], search_term)
                    if res_lhs isa Bool
                        res_rhs = check_rhs(eqns[i], search_term)
                        if !(res_rhs isa Bool)
                            sub_dict[res_rhs.lhs] = res_rhs.rhs
                            popat!(eqns, i)
                            break
                        end
                    else 
                        sub_dict[res_lhs.lhs] = res_lhs.rhs
                        popat!(eqns, i)
                        break
                    end
                end
            end
            eqn = simplify(substitute(eqn, sub_dict))
            eqns = map(x -> expand_derivatives(expand(simplify(substitute(x, sub_dict), expand = true))), eqns)
            @show eqn.rhs
            args = get_args(eqn.rhs)
            final_args = []
            for arg ∈ args
                check = 0
                for sv ∈ state_vars
                    if !((arg - sv == 0) isa Bool)
                        check += 1
                    end
                end
                (check == length(state_vars)) ? push!(final_args, arg) : ()
            end
        search_terms = union(search_terms, final_args)
        end
        eqns = map(x -> expand_derivatives(expand(simplify(substitute(x, Dict([eqn.lhs => eqn.rhs])),  expand = true))), eqns)
        eqn = expand(simplify(substitute(eqn, sub_dict)))
        push!(res_eqns, eqn)
    end
    return res_eqns
end

function resolve_derivative_causality!(BG::BondGraph)
    # Find the Elements with Derivative Causality aka element.causality=true
    # constiuitive_equations = []
    eqns = equations(BG.model)
    @show length(eqns)
    # Get q-variables
    fn(g, v) = (get_prop(g, v, :type) == :C || get_prop(g, v, :type) == :I) && !get_prop(g, v, :causality)
    state_vars = map(v -> get_prop(BG.graph, v, :state_var), filter_vertices(BG.graph, fn))
    state_vars = reduce(vcat, state_vars)
    @show state_vars
    # state_vars = reduce(vcat, map(x -> x.state_var, filter(x -> (x.causality == false), collect(values(BG.elements)))))
    D = Differential(BG.model.iv)
    diff_indexes = []
    implicit_eqns = []
    causal_nodes = filter_vertices(BG.graph, :causality, :true)
    # causal_names = map(v->get_prop(BG, v, :name), causal_nodes)
    for node ∈ causal_nodes
        name = get_prop(BG.graph, node, :name)
        @show name
        if get_prop(BG.graph, node, :type) == :I
            constiuitive_equation = BG[name].p ~ BG[name].I * BG[name].f
            old_eqn = BG[name].f ~ BG[name].p / BG[name].I
            index_old_eqn = indexin([old_eqn], eqns)[1]
            diff_eqns, alg_eqns = get_diff([eqns[1:(index_old_eqn - 1)];eqns[(index_old_eqn + 1):end]])
            @show constiuitive_equation
            imp_eqn = get_implicit(state_vars, [constiuitive_equation], alg_eqns)
            e_eqn = BG[name].e ~ expand(expand_derivatives(D(imp_eqn[1].rhs)))
            e_eqn = substitute(e_eqn, Dict(map(j -> diff_eqns[j].lhs => diff_eqns[j].rhs, eachindex(diff_eqns))))
            eqns[index_old_eqn] = old_eqn
            diff_eqn_index = indexin([D(BG[name].p) ~ BG[name].e], eqns)[1]
            push!(diff_indexes, diff_eqn_index)
            push!(implicit_eqns, e_eqn)
        elseif get_prop(BG.graph, node, :type) == :C
            constiuitive_equation = BG[name].q ~ BG[name].C * BG[name].e
            old_eqn = BG[name].e ~ BG[name].q / BG[name].C
            index_old_eqn = indexin([old_eqn], eqns)[1]
            diff_eqns, alg_eqns = get_diff([eqns[1:(index_old_eqn - 1)];eqns[(index_old_eqn + 1):end]])
            imp_eqn = get_implicit(state_vars, [constiuitive_equation], alg_eqns)
            f_eqn = BG[name].f ~ expand(expand_derivatives(D(imp_eqn[1].rhs)))
            f_eqn = substitute(f_eqn, Dict(map(j -> diff_eqns[j].lhs => diff_eqns[j].rhs, eachindex(diff_eqns))))
            eqns[index_old_eqn] = old_eqn
            diff_eqn_index = indexin([D(BG[name].q) ~ BG[name].f], eqns)[1]
            push!(diff_indexes, diff_eqn_index)
            push!(implicit_eqns, f_eqn)
        end
    end
    map(i -> eqns[diff_indexes[i]] = implicit_eqns[i], eachindex(implicit_eqns))
    BG.model = ODESystem(eqns, BG.model.iv, states(BG.model), parameters(BG.model))
    nothing
end