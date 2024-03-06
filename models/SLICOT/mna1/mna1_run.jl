using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@current_path "mna1_model.jl")
include(@current_path "mna1_specifications.jl")

S = mna1_model()
X0, options = mna1_specification()
options = Dict(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [i:i for i in 1:578]

# reachability algorithms
algorithm_reach = BFFPSV18(:δ => 1e-3, :vars => [1], :partition => partition)
algorithm_check = BFFPSV18(:δ => 4e-4, :vars => [1], :partition => partition)

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
savefig("mna1")
