using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "pde_model.jl")
include(@relpath "pde_specifications.jl")

S = pde_model()
X0, options = pde_specification()

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [i:i for i in 1:84]

# reachability algorithms
algorithm_reach = BFFPSV18(:δ => 1e-3, :vars => [1], :partition => partition)
algorithm_check = BFFPSV18(:δ => 3e-4, :vars => 1:84, :partition => partition)

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
savefig("pde")