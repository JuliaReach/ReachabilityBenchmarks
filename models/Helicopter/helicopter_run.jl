using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@current_path "helicopter_model.jl")
include(@current_path "helicopter_specifications.jl")

S = helicopter_model()
X0, options = helicopter_specification()
options = Options(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [[i] for i in 1:28]

# reachability algorithms
algorithm_reach = BFFPSV18(:δ => 5e-3, :vars => [1], :partition => partition,
                           :assume_homogeneous => true)
algorithm_check = BFFPSV18(:δ => 5e-3, :vars => [1], :partition => partition,
                           :assume_homogeneous => true)

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
savefig("helicopter")
