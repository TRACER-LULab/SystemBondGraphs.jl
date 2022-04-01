#=
function IP(;Dict(node=>variable))
    some math occurs here
    [out1, out2, out3]
end

=#
function add_IP!(BG::AbstractBondGraph, inputs, outputs, IP, name::Symbol)
    # Add IP to Graph
    add_vertex!(BG.graph)
    node_index = nv(BG.graph)
    set_prop!(BG.graph, node_index, :name, name)
    for (input, _) in inputs
        add_edge!(BG.graph, BG.graph[input, :name], BG.graph[name, :name])
    end
    for (output, _) in outputs
        add_edge!(BG.graph, BG.graph[name, :name], BG.graph[output, :name])
    end
    # Create outputs for information processor
    # inputs = ParentScope.(|)
    for (k, v) in inputs
        inputs[k] = ParentScope.(v)
    end
    # for (k, v) in outputs
    #     display(v)
    #     outputs[k] = ParentScope.(v)
    # end
    out_results = IP(inputs)
    subs = Pair{Any, Any}[]
    for (k, v) in out_results
        push!(subs,(outputs[k] .=> v)...)
    end

    # sys = ODESystem(eqns, BG.model.iv, [], [], name = name)
    props = Dict(
        :type => :IP, 
        :subs => subs
    )
    set_props!(BG.graph, node_index, props)
    subs
    # nothing
end