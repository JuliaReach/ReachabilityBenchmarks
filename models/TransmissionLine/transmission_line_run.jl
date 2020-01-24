using ReachabilityBenchmarks, MathematicalSystems, Reachability, Plots

include(@relpath "transmission_line_model.jl")
include(@relpath "transmission_line_specifications.jl")

# number of circuit nodes
η = 2

S = transmission_line_model(η)
X0, options = transmission_line_specification(S)

# initial value problem
problem = InitialValueProblem(S, X0)

# partition with 1D blocks
partition = [i:i for i in 1:(2 * η)]
partition = [1:(2 * η)]  # TODO temporarily use most precise partition

# variables of interest: U₁ ~ 1, U_η ~ η, I₁ ~ η+1, I_η ~ 2*η
vars = [1, η, η + 1, 2 * η]

# reachability algorithm
algorithm_reach = ASB07_decomposed(:δ => 0.002, :max_order => 400,
    :order_discretization => 9, :partition => partition, :vars => vars)
# TODO ":order_discretization => 6" was used in the paper but that did not work
# TODO temporary alternative with ASB07
algorithm_reach2 = ASB07(:δ => 0.002, :max_order => 400,
    :order_discretization => 9)
# TODO temporary alternative with BFFPSV18 without intervals
problem3 = InitialValueProblem(ConstrainedLinearControlContinuousSystem(
    mid(S.A), mid(S.B), nothing, S.U), X0)
algorithm_reach3 = BFFPSV18(:δ => 0.002, :partition => partition, :vars => vars)

# compute flowpipe
options[:mode] = "reach"
solution = solve(problem, options; op=algorithm_reach)
solution2 = solve(problem, options; op=algorithm_reach2)
solution3 = solve(problem3, options; op=algorithm_reach3)

# TODO temporary filtering of flowpipe to the first k sets
k = 10
flowpipes = [Flowpipe(solution.flowpipes[1].reachsets[1:k])]
solution = ReachSolution(flowpipes, solution.options)
flowpipes = [Flowpipe(solution2.flowpipes[1].reachsets[1:k])]
solution2 = ReachSolution(flowpipes, solution2.options)
flowpipes = [Flowpipe(solution3.flowpipes[1].reachsets[1:k])]
solution3 = ReachSolution(flowpipes, solution3.options)

for (vars, suffix) in [
                       ([0, η], "t-U_n"),
                       ([1, η + 1], "U_1-I_1"),
                       ([η, 2 * η], "U_n-I_n")
                      ]
    # project flowpipe to time and x_η
    solution.options[:plot_vars] = vars
    solution_proj = project(solution)
    solution2_proj = project(solution2)
    solution3_proj = project(solution3)

    # plot projection
    plot(solution_proj, opacity=0.4, color=:blue)
    plot!(solution2_proj, opacity=0.4, color=:red)
    plot!(solution3_proj, opacity=0.1, color=:green)
    savefig("transmission_line_" * suffix)

    # TODO temporarily also plot X0
    using SparseArrays
    if vars[1] == 0
        M = sparse([1], [vars[2]], [1.], 1, 2 * η)
        Y = LazySets.Interval(0., algorithm_reach.options[:δ]) × (M * X0)
    else
        M = sparse([1, 2], vars, [1., 1.], 2, 2 * η)
        Y = M * X0
    end
    plot!(Y, opacity=0.8, color=:yellow)
    savefig("transmission_line_X0_" * suffix)
end
