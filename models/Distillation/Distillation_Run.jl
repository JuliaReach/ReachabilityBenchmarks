using Reachability

problem = include("Distillation.jl")

# general options
𝑂_Dist = Options(:T => 6.0, :mode => "reach", :plot_vars => [1, 2])

# algorithm specific options
𝑂_Dist_BFFPSV18 = Options(:vars => [1, 2], :partition => [1:2, 3:4, 5:6, 7:8],
                          :δ => 0.003)

sol = solve(problem, 𝑂_Dist, op=BFFPSV18(𝑂_Dist_BFFPSV18))
