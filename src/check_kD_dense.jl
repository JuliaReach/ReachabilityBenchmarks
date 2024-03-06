# check property, kD blocks, intervals, given delta, dense time
using Reachability

check_kD_dense(model::String, delta::Float64) = check_kD_dense([model], [delta])

function check_kD_dense(models::Vector{String}, deltas::Vector{Float64})
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

    println("-- benchmark suite 'check_kD_dense' --")

    @assert length(models) == length(deltas) "the input lengths must be equal"
    for i in eachindex(models)
        model = models[i]
        delta = deltas[i]

        println("- analyzing model '$model' -")

        # raw model
        func = getfield(Main, Symbol(model))
        S, options_raw = func(:mode => "check")

        # options
        dict_raw = options_raw.dict
        dict_raw[:verbosity] = "info"
        dict_raw[:δ] = delta

        for i in 1:2
            dict = copy(dict_raw)
            if i == 1
                # warm-up run
                dict[:N] = 3
            else
                # benchmark settings
                dict[:T] = 20.0
                dict[:logfile] = "$model-check-kD-dense-givenstep.txt"
            end
            solve(S, Options(dict))
        end
        println()
    end
end
