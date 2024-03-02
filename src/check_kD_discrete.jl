# check property, kD blocks, intervals, fixed delta, discrete time
using Reachability

check_kD_discrete(model::String) = check_kD_discrete([model])

function check_kD_discrete(models::Vector{String})
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

    println("-- benchmark suite 'check_kD_discrete' --")

    for model in models
        println("- analyzing model '$model' -")

        # raw model
        func = getfield(Main, Symbol(model))
        S, options_raw = func(:mode => "check")

        # options
        dict_raw = options_raw.dict
        dict_raw[:verbosity] = "info"
        dict_raw[:δ] = 5e-3
        dict_raw[:approx_model] = "nobloating"
        dict_raw[:eager_checking] = false

        for i in 1:2
            dict = copy(dict_raw)
            if i == 1
                # warm-up run
                dict[:N] = 3
            else
                # benchmark settings
                dict[:T] = 20.0
                dict[:logfile] = "$model-check-kD-discrete-fixedstep.txt"
            end
            solve(S, Options(dict))
        end
        println()
    end
end
