# reach, 1D blocks, intervals, fixed delta, single variable
using Reachability

reach_1D_single(model::String) = reach_1D_single([model])

function reach_1D_single(models::Vector{String})
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
    if models_loaded
        println("models have been loaded; please run again for analysis")
        return
    end

    println("-- benchmark suite 'reach_1D_single' --")

    for model in models
        println("- analyzing model '$model' -")

        # raw model
        func = getfield(Main, Symbol(model))
        S, options_raw = func(:mode => "reach")

        # options
        n = Systems.statedim(S)
        dict_raw = options_raw.dict
        dict_raw[:verbosity] = "info"
        dict_raw[:δ] = 1e-3
        dict_raw[:partition] = [[i] for i in 1:n]
        dict_raw[:set_type] = Interval

        for i in 1:2
            dict = copy(dict_raw)
            if i == 1
                # warm-up run
                dict[:N] = 3
            else
                # benchmark settings
                dict[:T] = 20.
                dict[:logfile] = "$model-reach-1D-fixedstep-onevar.txt"
            end
            result = solve(S, Options(dict))
        end
        println()
    end
end
