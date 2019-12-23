using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "iss_model.jl")
include(@relpath "iss_specifications.jl")

S = iss_model()
X0, options = iss_specification()

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with mixed 1D and 135D blocks
partition = vcat([[i] for i in 1:135], [136:270])

# reachability algorithms
algorithm_reach = BFFPSV18(:δ => 1e-3, :vars => [182], :partition => partition,
                           :assume_sparse => true)
algorithm_check = BFFPSV18(:δ => 6e-4, :vars => 136:270, :partition => partition,
                           :assume_sparse => true, :lazy_inputs_interval => -1)

# check property
options[:mode] = "check"
solution = solve(problem, options; op=algorithm_check)
@assert solution.satisfied

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options; op=algorithm_reach)

# project flowpipe to time and x₁₈₂
solution.options[:plot_vars] = [0, 182]
solution_proj = project(solution)

# plot projection
plot(solution_proj)
savefig("iss")
