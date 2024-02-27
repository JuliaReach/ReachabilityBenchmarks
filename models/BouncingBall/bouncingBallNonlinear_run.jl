using ReachabilityBenchmarks, MathematicalSystems, Reachability,
      TaylorIntegration, Plots

include(@relpath "bouncingBallNonlinear_model.jl")
include(@relpath "bouncingBallNonlinear_specifications.jl")

S = bouncingBallNonlinear_model()
X0, options = bouncingBallNonlinear_specification()
options = Options(options)

# initial value problem
problem = InitialValueProblem(S, X0)

# reachability algorithm
algorithm_continuous = TMJets(:orderT => 5, :orderQ => 2, :abs_tol => 1e-10)
algorithm_hybrid = LazyDiscretePost(:check_invariant_intersection => true)

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options, algorithm_continuous, algorithm_hybrid)

# project flowpipe to x₁ and x₂
solution.options[:plot_vars] = [1, 2]
solution_proj = project(solution)

# plot projection
plot(solution_proj; xlab="x", ylab="v")
savefig("bouncingBallNonlinear")
