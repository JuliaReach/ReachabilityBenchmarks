using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "crane_model.jl")
include(@relpath "crane_specifications.jl")

S = crane_model()
X0, options = crane_specification()
options = Options(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [[i] for i in 1:6]

# reachability algorithms
algorithm_reach = BFFPSV18(:δ => 0.01, :vars => [1], :partition => partition)
algorithm_check = BFFPSV18(:δ => 0.01, :vars => [1], :partition => partition)

# check property
options[:mode] = "check"
solution = solve(problem, options; op=algorithm_check)
@assert solution.satisfied

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options; op=algorithm_reach)

# project flowpipe to time and x₁
solution.options[:plot_vars] = [0, 1]
solution_proj = project(solution)

# plot projection
plot(solution_proj)
savefig("crane")
