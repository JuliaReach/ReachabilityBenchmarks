using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@current_path "heat_model.jl")
include(@current_path "heat_specifications.jl")

S = heat_model()
X0, options = heat_specification()
options = Dict(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [i:i for i in 1:200]

# reachability algorithms
algorithm_reach = BFFPSV18(:δ => 1e-3, :vars => [133], :partition => partition)
algorithm_check = BFFPSV18(:δ => 1e-3, :vars => [133], :partition => partition)

# check property
options[:mode] = "check"
solution = solve(problem, options; op=algorithm_check)
@assert solution.satisfied

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options; op=algorithm_reach)

# project flowpipe to time and x₁₃₃
solution.options[:plot_vars] = [0, 133]
solution_proj = project(solution)

# plot projection
plot(solution_proj)
savefig("heat")
