using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "fom_model.jl")
include(@relpath "fom_specifications.jl")

S = fom_model()
X0, options = fom_specification()
options = Dict(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [i:i for i in 1:1006]

# reachability algorithms
algorithm_reach = BFFPSV18(:δ => 1e-3, :vars => [1], :partition => partition,
                           :exp_method => "lazy", :assume_sparse => true)
algorithm_check = BFFPSV18(:δ => 1e-3, :vars => 1:1006, :partition => partition,
                           :exp_method => "lazy", :assume_sparse => true,
                           :lazy_inputs_interval => -1)

# check property
# only verify a small time horizon (full time horizon does not verify)
options_check = copy(options)
options_check[:mode] = "check"
options_check[:T] = 0.001
solution = solve(problem, options_check; op=algorithm_check)
@assert solution.satisfied

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options; op=algorithm_reach)

# project flowpipe to time and x₁
solution.options[:plot_vars] = [0, 1]
solution_proj = project(solution)

# plot projection
plot(solution_proj)
savefig("fom")
