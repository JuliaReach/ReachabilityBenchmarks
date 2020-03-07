using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "fiveDimSys_model.jl")
include(@relpath "fiveDimSys_specifications.jl")

S = fiveDimSys_model()
X0, options = fiveDimSys_specification()
options = Options(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [[i] for i in 1:5]

# reachability algorithms
algorithm_reach1 = BFFPSV18(:δ => 0.005, :vars => [1, 2],
                            :partition => partition)
algorithm_reach2 = GLGM06(:δ => 0.005, :max_order => 5)

# compute flowpipe
options[:mode] = "reach"
solution1 = solve(problem, options; op=algorithm_reach1)
solution2 = solve(problem, options; op=algorithm_reach2)

# project flowpipe to time and x₁
solution1.options[:plot_vars] = [1, 2]
solution_proj1 = project(solution1)
solution2.options[:plot_vars] = [1, 2]
solution_proj2 = project(solution2)

# plot projection
plot(solution_proj1)
savefig("fiveDimSys1")
plot(solution_proj2)
savefig("fiveDimSys2")
