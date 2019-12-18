using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "beam_model.jl")
include(@relpath "beam_specifications.jl")

S = beam_model()
X0, options = beam_specification()

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [i:i for i in 1:348]

# reachability algorithms
algorithm_reach = BFFPSV18(:δ => 1e-3, :vars => [89], :partition => partition)
algorithm_check = BFFPSV18(:δ => 5e-5, :vars => [89], :partition => partition)

# check property
options[:mode] = "check"
solution = solve(problem, options; op=algorithm_check)
@assert solution.satisfied

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options; op=algorithm_reach)

# project flowpipe to time and x₈₉
solution.options[:plot_vars] = [0, 89]
solution_proj = project(solution)

# plot projection
plot(solution_proj)
savefig("beam")
