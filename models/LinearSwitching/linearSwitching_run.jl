using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "linearSwitching_model.jl")
include(@relpath "linearSwitching_specifications.jl")

S = linearSwitching_model()
X0, options = linearSwitching_specification()
options = Options(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# reachability algorithm
algorithm_continuous = BFFPS19(:δ => 1e-4, :partition => [1:2, 3:3, 4:4, 5:5])
algorithm_hybrid = DecomposedDiscretePost(:out_vars => [1, 2], :clustering => :none)

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options, algorithm_continuous, algorithm_hybrid)

# project flowpipe to x₁ and x₂
solution.options[:plot_vars] = [1, 2]
solution_proj = project(solution)

# plot projection
plot(solution_proj)
savefig("linearSwitching")
