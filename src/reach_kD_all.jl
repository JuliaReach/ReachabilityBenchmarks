# reach, varying block size, fixed delta, all variables
using Reachability

reach_kD_all(model::String, create_plots::Bool=false) =
    reach_kD_all([model], create_plots)

function reach_kD_all(models::Vector{String}, create_plots::Bool=false)
    # load models
    models_loaded = false
    for model in models
        if !isdefined(Main, Symbol(model))
            if !models_loaded
                model_library = include(@relpath "model_library.jl")
            end
            path = (@relpath "../") * model_library[model]
            include(path)
            models_loaded = true
        end
    end
    plots_loaded = false
    if create_plots && !isdefined(Main, :plot_reach)
        include(@relpath "plot_reach.jl")
        plots_loaded = true
    end
    if models_loaded
        println("models have been loaded; please run again for analysis")
        return
    elseif plots_loaded
        println("plot backend has been loaded; please run again for analysis")
        return
    end

    println("-- benchmark suite 'reach_kD_all' --")

    for model in models
        println("- analyzing model '$model' -")

        # raw model (check mode to obtain the property)
        func = getfield(Main, Symbol(model))
        S, options_raw = func(:mode => "check")
        n = MathematicalSystems.statedim(S)

        # extract output function from property
        property = options_raw[:property]
        if !(property isa LinearConstraintProperty)
            println("cannot construct output function from properties of " *
                "type $(typeof(property))... skipping model")
            break
        end
        projection_matrix = reshape(property.clauses[1].atoms[1].a, (1, n))

        # options
        dict_raw = options_raw.dict
        dict_raw[:mode] = "reach"
        dict_raw[:verbosity] = "info"
        dict_raw[:δ] = 1e-3
        dict_raw[:vars] = 1:n
        dict_raw[:plot_vars] = [0, n]
        dict_raw[:lazy_X0] = true
        dict_raw[:lazy_inputs_interval] = -1
        dict_raw[:project_reachset] = true
        dict_raw[:projection_matrix] = projection_matrix

        # create uniform partitions
        partitions = Vector{Vector{AbstractVector{Int}}}(ceil(Int, log(2, n)) + 1)
        k = 0
        while true
            m = 2^k
            k += 1
            if m >= n
                m = n
            end
            uniform = [m*i-(m-1):m*i for i in 1:div(n, m)]
            partitions[k] = n%m == 0 ? uniform : vcat(uniform, [div(n, m) *m+1:n])
            if m == n
                break
            end
        end

        for i in 1:(length(partitions) + 1)
            dict = copy(dict_raw)
            if i == 1
                # warm-up run
                dict[:N] = 3
                dict_raw[:partition] = partitions[1]
            else
                # benchmark settings
                k = i == length(partitions) + 1 ? n : 2^(i-2)
                dict[:N] = 50
                dict[:logfile] = "$model-reach-$(k)D-varying-fixedstep-allvars.txt"
                dict_raw[:partition] = partitions[i-1]
            end
            result = solve(S, Options(dict))
            if create_plots && i > 1
                plot_reach(result, "$model-$(k)D-varying-all")
            end
        end
        println()
    end
end
