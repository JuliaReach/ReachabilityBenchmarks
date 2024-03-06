# check property, 1D blocks, intervals, fixed delta, discrete time
using Reachability

check_1D_discrete(model::String) = check_1D_discrete([model])

function check_1D_discrete(models::Vector{String})
    # load models
    models_loaded = false
    for model in models
        if !isdefined(Main, Symbol(model))
            if !models_loaded
                model_library = include(@current_path "model_library.jl")
            end
            path = (@current_path "../") * model_library[model]
            include(path)
            models_loaded = true
        end
    end
    if models_loaded
        println("models have been loaded; please run again for analysis")
        return
    end

    println("-- benchmark suite 'check_1D_discrete' --")

    for model in models
        println("- analyzing model '$model' -")

        # raw model
        func = getfield(Main, Symbol(model))
        S, options_raw = func(:mode => "check")

        # options
        n = MathematicalSystems.statedim(S)
        dict_raw = options_raw.dict
        dict_raw[:verbosity] = "info"
        dict_raw[:δ] = 5e-3
        dict_raw[:approx_model] = "nobloating"
        dict_raw[:partition] = [[i] for i in 1:n]
        dict_raw[:set_type] = Interval
        dict_raw[:lazy_inputs_interval] = 0
        dict_raw[:eager_checking] = false

        for i in 1:2
            dict = copy(dict_raw)
            if i == 1
                # warm-up run
                dict[:N] = 3
            else
                # benchmark settings
                dict[:T] = 20.0
                dict[:logfile] = "$model-check-1D-discrete-fixedstep.txt"
            end
            solve(S, Options(dict))
        end
        println()
    end
end
