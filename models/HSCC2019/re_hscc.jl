using HybridSystems, MathematicalSystems, LazySets, Reachability, Polyhedra, Optim

import LazySets.HalfSpace
import LazySets.Approximations: overapproximate, OctDirections

include("FilteredOscillator.jl")

function get_projection(sol::AbstractSolution, system_dimension::Int64)::AbstractSolution
    N = Float64
    sol_processed =  Reachability.ReachSolution(
        [Reachability.ReachSet{CartesianProductArray{N}, N}(
            CartesianProductArray{N, HPolytope{N}}(
                [overapproximate(rs.X, OctDirections(system_dimension))]),
            rs.t_start, rs.t_end) for rs in sol.Xk],
        sol.options)

    sol_proj = Reachability.ReachSolution(Reachability.project_reach(
        sol_processed.Xk, [1,2], system_dimension, sol.options), sol.options);

    return sol_proj;
end

function warmup(opDs)
    # warmup run for each opD for low dimension
    println("warm-up runs")
    run (2, opDs, false, nothing)
    run (2, opDs, false, nothing)
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
            sol = filtered_oscillator(n0, opD, T);
        end
        if project_and_store
            sol_proj = get_projection(sol, n0+3);
            push!(results, (sol_proj, name));
        end
    end
end

function benchmark(project_and_store::Bool=false)
    # discrete-post operators + short name + upper bound on dimensionality
    opDs = [
            (ConcreteDiscretePost(),      "C",   8)
            (LazyDiscretePost(),          "L", 8)
            (ApproximatingDiscretePost(), "A", 8)
           ];

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
            n0 *= 2;
        end
    end

    return results
end

return benchmark(true)
