using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "transmission_line_model.jl")
include(@relpath "transmission_line_specifications.jl")

# number of circuit nodes
η = 20

S = transmission_line_model(η)
X0, options = transmission_line_specification(S)

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [i:i for i in 1:(2 * η)]

# variables of interest: U_η ~ η, I_η ~ 2*η, U₁ ~ 1, I₁ ~ η+1, 
vars = [1, η, η + 1, 2 * η]

# reachability algorithm
algorithm_reach = ASB07(:δ => 0.002, :max_order => 400,
                        :order_discretization => 6, :partition => partition,
                        :vars => vars)

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options; op=algorithm_reach)

for (vars, suffix) in [([0, η], "t-U_n"),
                       ([0, η], "U_1-I_n"),
                       ([0, η], "U_n-I_n")]
    # project flowpipe to time and x_η
    solution.options[:plot_vars] = vars
    solution_proj = project(solution)

    # plot projection
    plot(solution_proj)
    savefig("transmission_line_" * suffix)
end
