using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "lorenz_model.jl")
include(@relpath "lorenz_specifications.jl")

S = lorenz_model()
X0, options = lorenz_specification()
options = Options(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# reachability algorithm
algorithm_reach = TMJets(:abs_tol=>1e-27, :orderT=>10, :orderQ=>2,
                         :max_steps=>50_000)

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options; op=algorithm_reach)

# project flowpipe to x₁ and x₃
solution.options[:plot_vars] = [1, 3]
solution_proj = project(solution)

# plot projection
plot(solution_proj)
savefig("lorenz")
