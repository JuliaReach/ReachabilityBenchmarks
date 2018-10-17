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
opDs = [Reachability.ReachSets.TextbookDiscretePost(),
            Reachability.ReachSets.LazyTextbookDiscretePost(),
            Reachability.ReachSets.ApproximatingDiscretePost()];

#warmup run for each opD for low dimension
println("Warmup run just after restart REPL")
for opD_i in 1:length(opDs)
    n0 = 2;
    while (n0 <= 4)
        opD = opDs[opD_i]
        sol = filtered_oscillator(n0, opDs[opD_i], 20., 20);
        n0 = n0*2;
    end
end
println("End of warmup run")

results = Vector{Tuple{AbstractSolution, Int64}}()
for opD_i in 1:length(opDs)
    println("**********************************")
    println(opDs[opD_i])
    upper_bound = opD_i == 1 ? 4 : 256;
    n0 = 2;
    while (n0 <= upper_bound)
        if n0 == 2
            t = 20.
            max_jumps = 20
        elseif n0 == 4
            t = 20.
            max_jumps = 20
        elseif n0 == 8
            t = 99.
            max_jumps = 20
        elseif n0 == 16
            t = 99.
            max_jumps = 30
        elseif n0 == 32
            t = 99.
            max_jumps = 40
        elseif n0 == 64
            t = 99.
            max_jumps = 40
        elseif n0 == 128
            t = 99.
            max_jumps = 1000
        elseif n0 == 196
            t = 99.
            max_jumps = 1000
        elseif n0 == 256
            t = 99.
            max_jumps = 1000
        end
        println("\t", n0)
        opD = opDs[opD_i];
        sol = filtered_oscillator(n0, opD, t, max_jumps);
        sol_proj = get_projection(sol, n0+2);
        push!(results, (sol_proj, opD_i));
        if n0 == 128
            n0 = 196
        elseif n0 == 196
            n0 = 256
        else
            n0 *= 2;
        end
    end
end

return results
