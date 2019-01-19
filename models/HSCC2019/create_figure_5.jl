# To recreate Figure 3, run the following code from the REPL.
# Just type:
#
#     include("create_figure_5.jl")

using HybridSystems, MathematicalSystems, LazySets, Reachability, Polyhedra, Optim, Plots
import LazySets.HalfSpace

#disable graphics output in GR - https://github.com/JuliaPlots/Plots.jl/issues/1182
ENV["GKSwstype"] = "100"

include("FilteredOscillator.jl")

N = Float64

function get_projection(sol::AbstractSolution, system_dimension::Int64) ::AbstractSolution
    sol_processed =  Reachability.ReachSolution([Reachability.ReachSet{CartesianProductArray{N}, N}(
                CartesianProductArray{N, HPolytope{N}}([LazySets.Approximations.overapproximate(rs.X, LazySets.Approximations.OctDirections(system_dimension))]),
                rs.t_start, rs.t_end) for rs in sol.Xk], sol.options)

    sol_proj = Reachability.ReachSolution(Reachability.project_reach(
        sol_processed.Xk, [1,system_dimension - 1], system_dimension, sol.options), sol.options);

    return sol_proj;
end

sol4 = filtered_oscillator(4, Reachability.ReachSets.ApproximatingDiscretePost(), 20., 20);
sol_proj4 = get_projection(sol4, 7);
plot(sol_proj4, tickfont = font(14))
savefig("Figure 5 Left.pdf")

sol16 = filtered_oscillator(16, Reachability.ReachSets.ApproximatingDiscretePost(), 99., 30);
sol_proj16 = get_projection(sol16, 19);
plot(sol_proj16, tickfont = font(14))
savefig("Figure 5 Right.pdf")
