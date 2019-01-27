include("FilteredOscillator.jl")

function warmup(opDs)
    println("warm-up runs")
    run(2, opDs, false, nothing)
    run(2, opDs, false, nothing)
    println("end of warm-up runs")
end

function run(n0, opDs, project_and_store, results)
    for (opD, name, upper_bound) in opDs
        if n0 > upper_bound
            continue
        end
        println("**********************************")
        println("$name $(n0)")
        if n0 <= 4
            T = 20.
        else
            T = 99.
        end
        @time begin
            sol = filtered_oscillator(n0, opD, T)
        end
        if project_and_store
            sol_proj = get_projection(sol, n0+3, [1, 2])
            push!(results, (sol_proj, name))
        end
    end
end

function benchmark(project_and_store::Bool=false)
    # discrete-post operators + short name + upper bound on dimensionality
    opDs = [
            (ConcreteDiscretePost(),      "C",   8)
            (LazyDiscretePost(),          "L", 256)
            (ApproximatingDiscretePost(), "A", 256)
           ]

    warmup(opDs)

    if project_and_store
        println("Note: projecting and returning output")
    else
        println("Note: skipping projection and storage output")
    end

    results = Vector{Tuple{AbstractSolution, String}}()
    n0 = 2
    while n0 <= 256
        run(n0, opDs, project_and_store, results)
        if n0 == 128
            n0 = 196
        elseif n0 == 196
            n0 = 256
        else
            n0 *= 2
        end
    end

    return results
end

return benchmark(false)
