function get_diff_and_alg_eqns(eqns::Vector{Equation})
    differential_equations = Equation[]
    algebraic_equations = Equation[]
    for eqn ∈ eqns
        if !(eqn.lhs isa Real) # Check if the LHS is not equal to zero
            if SymbolicUtils.operation(eqn.lhs) isa Differential
                push!(differential_equations, eqn)
            else
                push!(algebraic_equations, eqn)
            end
        else
            push!(algebraic_equations, eqn)
        end
    end
    return differential_equations, algebraic_equations
end

function terms_in_final(term_list, final_syms)
    complete = true
    for term ∈ term_list
        idx = indexin([term], final_syms)[1]
        complete = complete && (idx !== nothing)
    end
    return complete
end

function get_implicit_relationship(constituitive_equation::Equation, algebraic_equations::Vector{Equation}, differential_equations::Vector{Equation}, state_variables::Vector{Num}, params)
    substitution_dict = Dict([])
    old_length = -1
    copy_algebraic_equations = copy(algebraic_equations)
    check(x) = (isa(x, Term) && (indexin([x], params)[1] === nothing) && (indexin([x], state_variables)[1] === nothing))
    term_list = filter(check, get_variables(constituitive_equation.rhs))
    final_syms = [state_variables]
    term_list_old = []
    # Check if term_list is equivalent to final_syms
    while !terms_in_final(term_list, final_syms) && (length(copy_algebraic_equations) > 0)
        term_list_old = term_list
        old_length = length(keys(substitution_dict))
        for term ∈ term_list
            found = false
            for i ∈ eachindex(copy_algebraic_equations)
                # Check if term is in eqn
                eqn = copy_algebraic_equations[i]
                equation_variables = get_variables(eqn)
                # Check if term is in equation variables
                index = indexin([term], equation_variables)[1]
                # If term is in equation variables
                if index isa Real
                    # Sovle for term and update substitution dictionary
                    eqn = flatten_fractions(expand(eqn.rhs))
                    if hasfield(typeof(eqn), :num)
                        eqn = 0 ~ eqn.num
                    else
                        eqn = 0 ~ eqn
                    end
                    substitution_dict[term] = Symbolics.solve_for(eqn, term)
                    # Remove Equation from Search
                    popat!(copy_algebraic_equations, i)
                    # Exit for loop
                    found = true
                end
                (found) ? (break) : ()
            end
            (found) ? (break) : ()
        end
        # Update Constituitive Equations
        constituitive_equation = expand(simplify(substitute(constituitive_equation, substitution_dict)))
        # Update Algebraic Equations
        copy_algebraic_equations = map(x -> substitute(x, substitution_dict), copy_algebraic_equations)
        # Update Term list in Constituitive Equation
        term_list = filter(check, get_variables(constituitive_equation.rhs))
    end
    constituitive_equation
end

function remove_casuality(BG::AbstractBondGraph; skip::Vector{Num} = Num[])
    # Make Model from Graph
    model = BondGraphs.graph_to_model(BG)
    # model = extend(BG.model, model)
    # Get Equations
    eqns = full_equations(model)
    # Get parameters
    params = parameters(model)
    # Create Differential Operator
    D = Differential(BG.model.iv)
    # Storage Arrays
    diff_indexes = []
    implicit_eqns = []
    # Get State-Variables for Derivative Casuality Elements
    fn(g, v) = (get_prop(g, v, :type) == :C || get_prop(g, v, :type) == :I) && !get_prop(g, v, :causality)
    state_vars = map(v -> get_prop(BG.graph, v, :state_var), filter_vertices(BG.graph, fn))
    state_vars = reduce(vcat, state_vars)
    # Get Nodes with Derivative Causality
    dc_nodes = filter_vertices(BG.graph, :causality, :true)
    for node ∈ dc_nodes
        element_type = get_prop(BG.graph, node, :type)
        name = get_prop(BG.graph, node, :name)
        if element_type == :I
            consituitive_equation = BG[name].p ~ BG[name].I * BG[name].f
            starting_equation = BG[name].f ~ BG[name].p / BG[name].I
            start_index = indexin([starting_equation], eqns)[1]
            # Sort and rearrange equations minus the starting equation
            differential_eqns, algebraic_eqns = get_diff_and_alg_eqns([eqns[1:(start_index-1)]; eqns[(start_index+1):end]])
            implicit_equation = get_implicit_relationship(consituitive_equation, algebraic_eqns, differential_eqns, state_vars, params)
            e_eqn = BG[name].e ~ expand(expand_derivatives(D(implicit_equation.rhs)))
            e_eqn = substitute(implicit_equation,
                Dict(
                    map(
                        j -> differential_eqns[j].lhs => differential_eqns[j].rhs,
                        eachindex(differential_eqns)
                    )
                )
            )
            eqns[start_index] = starting_equation
            diff_eqn_index = indexin([D(BG[name].p) ~ BG[name].e], eqns)[1]
            push!(diff_indexes, diff_eqn_index)
            push!(implicit_eqns, e_eqn)
        elseif element_type == :C
            consituitive_equation = BG[name].q ~ BG[name].C * BG[name].e
            starting_equation = BG[name].e ~ BG[name].q / BG[name].C
            start_index = indexin([starting_equation], eqns)[1]
            # Sort and rearrange equations minus the starting equation
            differential_eqns, algebraic_eqns = get_diff_and_alg_eqns([eqns[1:(start_index-1)]; eqns[(start_index+1):end]])
            implicit_equation = get_implicit_relationship(consituitive_equation, algebraic_eqns, differential_eqns, state_vars, params)
            f_eqn = BG[name].f ~ expand(expand_derivatives(D(implicit_equation.rhs)))
            f_eqn = substitute(f_eqn,
                Dict(
                    map(
                        j -> differential_eqns[j].lhs => differential_eqns[j].rhs,
                        eachindex(differential_eqns)
                    )
                )
            )
            eqns[start_index] = starting_equation
            diff_eqn_index = indexin([D(BG[name].q) ~ BG[name].f], eqns)[1]
            push!(diff_indexes, diff_eqn_index)
            push!(implicit_eqns, f_eqn)
        end
    end
    map(i -> eqns[diff_indexes[i]] = implicit_eqns[i], eachindex(implicit_eqns))
    @named model = ODESystem(eqns)
end
