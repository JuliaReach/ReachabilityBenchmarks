using HybridSystems, MathematicalSystems, LazySets, Reachability, Polyhedra, Optim
import LazySets.HalfSpace
include("FilteredOscillator.jl")

N = Float64

function get_projection(sol::AbstractSolution, system_dimension::Int64) ::AbstractSolution
    sol_processed =  Reachability.ReachSolution([Reachability.ReachSet{CartesianProductArray{N}, N}(
                CartesianProductArray{N, HPolytope{N}}([LazySets.Approximations.overapproximate(rs.X, LazySets.Approximations.OctDirections(system_dimension))]),
                rs.t_start, rs.t_end) for rs in sol.Xk], sol.options)

    sol_proj = Reachability.ReachSolution(Reachability.project_reach(
        sol_processed.Xk, [1,2], system_dimension, sol.options), sol.options);

    return sol_proj;
end

# discrete post operators + short name + upper bound on dimensionality
opDs = [
        (Reachability.ReachSets.TextbookDiscretePost(),      "C",   8)
        (Reachability.ReachSets.LazyTextbookDiscretePost(),  "L", 256)
        (Reachability.ReachSets.ApproximatingDiscretePost(), "A", 256)
       ];

#warmup run for each opD for low dimension
println("Warmup run just after restart REPL")
for (opD, name, upper_bound) in opDs
    n0 = 2;
    while (n0 <= 4)
        sol = filtered_oscillator(n0, opD, 20., 20);
        n0 = n0*2;
    end
end
println("End of warmup run")

project_and_store = false
if project_and_store
    println("Note: using output projection and storage")
else
    println("Note: skipping output projection and storage")
end

results = Vector{Tuple{AbstractSolution, String}}()
n0 = 2
while n0 <= 256
    for (opD, name, upper_bound) in opDs
        println("**********************************")
        if n0 > upper_bound
            continue
        end
        println("$name $(n0)")
        if n0 == 2
            T = 20.
            max_jumps = 20
        elseif n0 == 4
            T = 20.
            max_jumps = 20
        elseif n0 == 8
            T = 99.
            max_jumps = 20
        elseif n0 == 16
            T = 99.
            max_jumps = 30
        elseif n0 == 32
            T = 99.
            max_jumps = 40
        elseif n0 == 64
            T = 99.
            max_jumps = 40
        elseif n0 == 128
            T = 99.
            max_jumps = 1000
        elseif n0 == 196
            T = 99.
            max_jumps = 1000
        elseif n0 == 256
            T = 99.
            max_jumps = 1000
        end
        sol = filtered_oscillator(n0, opD, T, max_jumps);
        if project_and_store
            sol_proj = get_projection(sol, n0+2);
            push!(results, (sol_proj, name));
        end
    end
    if n0 == 128
        n0 = 196
    elseif n0 == 196
        n0 = 256
    else
        n0 *= 2;
    end
end

return results
