"""

Remove Algebraic Constraint Equations Generated by ModelingToolkit

"""
# function remove_algebraic(BG::AbstractBondGraph, model)
#     ## Remove Algebraic Constraint Equation
#     # model = structural_simplify(model)
#     filter(eq -> eq.lhs isa Int64, full_equations(model))
#     eqns = full_equations(model)
#     subs_dict = Dict{Term, Any}()
#     obs_eqns = Equation[]
#     while length(filter(eq -> eq.lhs isa Int64, eqns)) > 0
#         state_var_nodes = filter_vertices(BG.graph, (g, v) -> has_prop(g, v, :state_var)) |> collect
#         state_vars = vcat(map(v -> get_prop(BG.graph, v, :state_var), state_var_nodes)...)
#         for eqn_index ∈ eachindex(eqns)
#             eqn = eqns[eqn_index]
#             if (eqn.lhs == 0) isa Bool
#                 vars = get_variables(eqn)
#                 terms = filter(x -> isa(x, Term), vars)
#                 # find the extra variables, not a state variable
#                 for term_index ∈ eachindex(terms)
#                     is_state = false
#                     for sv ∈ state_vars
#                         if isa(sv - terms[term_index] == 0, Bool)
#                             is_state = true
#                             break
#                         end
#                     end
#                     if !is_state
#                         eqn = flatten_fractions(expand(eqn.rhs))
#                         if hasfield(typeof(eqn), :den)
#                             eqn = 0 ~ eqn.num
#                         else
#                             eqn = 0 ~ eqn
#                         end
#                         res = Symbolics.solve_for(eqn, terms[term_index])
#                         subs_dict[terms[term_index]] = res
#                         # Apply to full_equations
#                         push!(obs_eqns, terms[term_index] ~ res)
#                         popat!(eqns, eqn_index)
#                         eqns = map(eqn -> substitute(eqn, subs_dict), eqns)
#                         eqns = simplify_fractions.(expand.(eqns))
#                         # Remove algebraic equation
#                         # Update subs_dict
#                         for (k, v) in pairs(subs_dict)
#                             subs_dict[k] = simplify_fractions(substitute(v, subs_dict))
#                         end
#                         break
#                     end
#                 end
#                 break
#             end
#         end
#     end
#     state_vars = vcat(get_variables.(getfield.(eqns, :lhs))...)
#     ode_model = ODESystem(eqns, model.iv, name = model.name)
#     return ode_model, obs_eqns
# end

"""
Traverse the BondGraph Structure to create an ODESystem
"""
function _generate_model(bg)
    # Check if graph is disconnected
    if !is_connected(bg.graph)
        error(
            "An element(s) in the bond graph is not connected to the rest of the graph. Please ensure that all elements are connected.",
        )
    end
    # Algebraic Relationships for the BondGraph
    eqns = Equation[]

    # Find all One Junction Nodes
    one_junctions = filter(x -> bg[x].type == :OneJunction, collect(labels(bg)))
    for J1 ∈ one_junctions
        out_nodes = outneighbor_labels(bg, J1)
        out_nodes = map(x -> (J1, x), out_nodes)
        # out_nodes = map(x ->x, out_nodes)

        in_nodes = inneighbor_labels(bg, J1) |> collect
        in_nodes = map(x -> (x, J1), in_nodes)
        # in_nodes = map(x->x, in_nodes)

        if !isempty(out_nodes)
            out_sum = sum(x -> bg[x...].model.e, out_nodes)
        else
            out_sum = 0
        end
        if !isempty(in_nodes)
            in_sum = sum(x -> bg[x...].model.e, in_nodes)
        else
            in_sum = 0
        end
        push!(eqns, 0 ~ in_sum - out_sum)
        nodes = [out_nodes; in_nodes]

        for i ∈ 2:length(nodes)
            push!(eqns, 0 ~ bg[nodes[i]...].model.f - bg[nodes[i-1]...].model.f)
            # push!(eqns, 0 ~ bg[nodes[i]].model.f - bg[nodes[i-1]].model.f)
        end
    end

    # Find all Zero Junction Nodes
    zero_junctions = filter(x -> bg[x].type == :ZeroJunction, collect(labels(bg)))

    for J0 ∈ zero_junctions
        out_nodes = outneighbor_labels(bg, J0)
        out_nodes = map(x -> (J0, x), out_nodes)

        in_nodes = inneighbor_labels(bg, J0) |> collect
        in_nodes = map(x -> (x, J0), in_nodes)
        if !isempty(out_nodes)
            out_sum = sum(x -> bg[x...].model.f, out_nodes)
        else
            out_sum = 0
        end
        if !isempty(in_nodes)
            in_sum = sum(x -> bg[x...].model.f, in_nodes)
        else
            in_sum = 0
        end
        push!(eqns, 0 ~ in_sum - out_sum)
        nodes = [out_nodes; in_nodes]
        nodes = map(x -> x, nodes)
        for i ∈ 2:length(nodes)
            push!(eqns, 0 ~ bg[nodes[i]...].model.e - bg[nodes[i-1]...].model.e)
        end
    end

    # Check for the two_port relationships
    two_ports = filter(x -> bg[x].type in [:Re, :TF, :GY, :MTF, :MGY], collect(labels(bg)))
    two_ports_sys = map(v -> bg[v].model, two_ports)

    for TP in two_ports
        ins = inneighbor_labels(bg, TP)
        for i in ins
            push!(eqns, 0 ~ bg[TP].model.f_in - bg[i, TP].model.f)
            push!(eqns, 0 ~ bg[TP].model.e_in - bg[i, TP].model.e)
        end
        outs = outneighbor_labels(bg, TP)
        for o in outs
            push!(eqns, 0 ~ bg[TP].model.f_out - bg[TP, o].model.f)
            push!(eqns, 0 ~ bg[TP].model.e_out - bg[TP, o].model.e)
        end
    end

    # Get all one-port systems
    elements = filter(
        x -> bg[x].type in [:B, :R, :C, :I, :M, :Ce, :Se, :Sf, :MPC, :MPI, :MPR],
        collect(labels(bg)),
    )
    for element in elements
        ins = inneighbor_labels(bg, element)
        for i in ins
            push!(eqns, 0 ~ bg[element].model.f - bg[i, element].model.f)
            push!(eqns, 0 ~ bg[element].model.e - bg[i, element].model.e)
        end
        outs = outneighbor_labels(bg, element)
        for o in outs
            push!(eqns, 0 ~ bg[element].model.f - bg[element, o].model.f)
            push!(eqns, 0 ~ bg[element].model.e - bg[element, o].model.e)
        end
    end
    element_sys = map(v -> bg[v].model, elements)
    @named system =
        ODESystem(eqns, bg.graph_data.iv, [], [], systems = [two_ports_sys; element_sys])
    return system
end

function generate_model(
    bg::MetaGraph{A,B,C,D,E,DATA,F,G},
) where {A,B,C,D,E,DATA<:BondGraphData,F,G}
    _generate_model(bg)
end
"""

Generate an ODE System from the BondGraph Structure

"""
function generate_model(
    bg::MetaGraph{A,B,C,D,E,DATA,F,G},
) where {A,B,C,D,E,DATA<:ChemBondGraphData,F,G}
    model = _generate_model(bg)
    model = structural_simplify(model)
    RW = SymbolicUtils.Rewriters
    r1 = @acrule log(~x) + log(~y) => log((~x) * (~y))
    r2 = @rule log(~x) - log(~y) => log((~x) / (~y))
    r3 = @rule (~x) * log(~y) => log((~y)^(~x))
    r4 = @rule exp(log(~x)) => ~x
    r5 = @acrule exp((~x) + (~y)) => exp(~x) * exp(~y)
    rw1 = RW.Fixpoint(RW.Chain([r1, r2, r3, r4, r5]))
    rw2 = RW.Prewalk(RW.Chain([r1, r2, r3, r4, r5]))
    rw3 = RW.Postwalk(RW.Chain([r1, r2, r3, r4, r5]))
    eqns = full_equations(model)
    for i ∈ eachindex(eqns)
        eqns[i] = eqns[i].lhs ~ eqns[i].rhs |> rw3 |> rw2 |> rw1 |> expand
    end
    defaults = model.defaults
    model = ODESystem(
        eqns,
        model.iv,
        unknowns(model),
        parameters(model),
        name = nameof(model),
        observed = observed(model),
        defaults = defaults,
    )
    return structural_simplify(model)
end

# function generate_model!(BG::BioBondGraph)
#     BG.model = generate_model(BG)
# end
