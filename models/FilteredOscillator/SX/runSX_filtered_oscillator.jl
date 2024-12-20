# Filtered oscillator
# ============================

using ReachabilityBenchmarks
using HybridSystems, MathematicalSystems, LazySets, Plots
import Polyhedra
using SpaceExParser: parse_sxmath, Basic, readsxmodel

AFFINE_SYSTEM = ConstrainedLinearControlContinuousSystem
HS = readsxmodel(@current_path "filtered_oscillator_flattened.xml"; ST=AFFINE_SYSTEM)

expr_p = parse_sxmath("0.2<=x & x<=0.3 & -0.1<=y & y<=0.1 & z==0 & x1==0 & x2==0 & x3==0")
#expr = [x,y,x1,x2,x3,z]
x0sets = []
vars = Basic[:x, :y, :x1, :x2, :x3, :z]
for expr_i in expr_p
    if LazySets.is_halfspace(expr_i)
        push!(x0sets, LazySets.convert(HalfSpace, expr_i; vars=vars))
    elseif LazySets.is_hyperplane(expr_i)
        push!(x0sets, LazySets.convert(Hyperplane, expr_i; vars=vars))
    end
end

prob = InitialValueProblem(HS, x0sets)

# calculate reachable states up to time T
input_options = Options(:mode => "reach");
problem_options = Options(:vars => 1:system_dimension, :T => 5.0, :δ => 0.05, :plot_vars => [1, 2],
                          :max_jumps => 1, :verbosity => 1, :init_loc => 3);
options_input = merge(problem_options, input_options)
sol = solve(prob, options_input);

plot(sol; indices=1:2:length(sol.Xk))
