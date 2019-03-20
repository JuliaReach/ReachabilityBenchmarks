using LazySets, Reachability, Polyhedra, Optim
import LazySets.HalfSpace
include("DriveTrain.jl");

# compatibility of julia versions
if VERSION >= v"0.7"
    using LinearAlgebra, SparseArrays
end

N = Float64

function get_projection(sol::AbstractSolution, system_dimension::Int64) ::AbstractSolution
    sol_processed =  Reachability.ReachSolution([Reachability.ReachSet{CartesianProductArray{N}, N}(
                CartesianProductArray{N, HPolytope{N}}([LazySets.Approximations.overapproximate(rs.X, LazySets.Approximations.OctDirections(system_dimension))]),
                rs.t_start, rs.t_end) for rs in sol.Xk], sol.options)

    sol_proj = Reachability.ReachSolution(Reachability.project_reach(
        sol_processed.Xk, [1,2], system_dimension, sol.options), sol.options);

    return sol_proj;
end

function get_solution(drivetrain, θ, opD, t, max_jumps)
    # initial condition in mode 1
    c = vcat([-0.0432, -11, 0, 30, 0, 30, 360], repeat([-0.0013, 30], θ))
    g = zeros(2*θ + 7,2*θ + 7)
    g[1,1], g[2,2], g[3,3], g[4,4], g[5,5], g[6,6], g[7,7] = 0.0056, 4.67, 0, 10, 0, 10, 120
    i = 1;
    while (i <= θ)
        g[2*i+6,2*i+6], g[2*i+7,2*i+7] = 0.0006, 10
        i+=1
    end
    X0 = Zonotope(c, g)

    system = InitialValueProblem(drivetrain, [(1, X0)]);
    plot_vars = [1, 3]
    system_dimension = 2*θ+7
    options = Options(:mode=>"reach",:vars=>1:system_dimension, :T=>t, :δ=>0.01,
                      :max_jumps=>max_jumps, :plot_vars=>plot_vars,
                      :ε_proj=>0.001, :verbosity=>0, :project_reachset=>false);


    # default algorithm
    @time begin
        sol = solve(system, options, Reachability.BFFPSV18(), opD);
    end
end

opDs = [
        (Reachability.ReachSets.LazyDiscretePost(),  "L", 64)
        (Reachability.ReachSets.ApproximatingDiscretePost(), "A", 64)
       ];


#warmup run for each opD for low dimension
println("Warmup run just after restart REPL")
for (opD, name, upper_bound) in opDs
    θ = 1;
    while (θ <= 2)
        HS = drivetrain(θ)
        sol = get_solution(HS, θ, opD, 2., 4);
        θ += 1;
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
θ = 1
while θ <= 47
    for (opD, name, upper_bound) in opDs
        println("**********************************")
        if θ > upper_bound
            continue
        end
        println("$name $(n0)")
        T = 2.
        max_jumps = 4
        end
        HS = drivetrain(θ)
        sol = get_solution(HS, θ, opD, T, max_jumps);
        if project_and_store
            sol_proj = get_projection(sol, 2*θ+7);
            push!(results, (sol_proj, name));
        end
    end
    θ += 1;
end

return results
