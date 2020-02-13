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
δ = 0.002
algorithm_reach = ASB07_decomposed(:δ => δ, :max_order => 400,
    :order_discretization => 9, :partition => partition, :vars => vars)
# TODO ":order_discretization => 6" was used in the paper but that did not work
# TODO temporary alternative with ASB07
algorithm_reach2 = ASB07(:δ => δ, :max_order => 400,
    :order_discretization => 9)
# TODO temporary alternative with BFFPSV18 without intervals
problem3 = InitialValueProblem(ConstrainedLinearControlContinuousSystem(
    mid(S.A), mid(S.B), nothing, S.U), X0)
algorithm_reach3 = BFFPSV18(:δ => δ, :partition => partition, :vars => vars)

# compute flowpipe
options[:mode] = "reach"
k = 50 # Int(div(options[:T], δ))
options[:T] = δ * k
solution = solve(problem, copy(options); op=algorithm_reach)
solution2 = solve(problem, copy(options); op=algorithm_reach2)
solution3 = solve(problem3, copy(options); op=algorithm_reach3)

for (vars, suffix) in [
                       ([0, η], "t-U_n"),
                       ([1, η + 1], "U_1-I_1"),
                       ([η, 2 * η], "U_n-I_n")
                      ]
    # project flowpipe to time and x_η
    solution.options[:plot_vars] = vars
    solution_proj = project(solution)
    solution2.options[:plot_vars] = vars
    solution2_proj = project(solution2)
    solution3.options[:plot_vars] = vars
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
        Y = LazySets.Interval(0., δ) × (M * X0)
    else
        M = sparse([1, 2], vars, [1., 1.], 2, 2 * η)
        Y = M * X0
    end
    plot!(Y, opacity=0.8, color=:yellow)
    savefig("transmission_line_X0_" * suffix)
end
