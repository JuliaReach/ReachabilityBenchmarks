using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "projectile_model.jl")
include(@relpath "projectile_specifications.jl")

S = projectile_model()
X0, options = projectile_specification()
options = Options(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 2D blocks
partition = [(2 * i - 1):(2 * i) for i in 1:2]

# reachability algorithms
algorithm_reach = BFFPSV18(:δ => 0.5, :vars => [1, 3], :partition => partition)

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options; op=algorithm_reach)

# project flowpipe to x₁ and x₃
solution.options[:plot_vars] = [1, 3]
solution_proj = project(solution)

# plot projection
plot(solution_proj)
savefig("projectile")
