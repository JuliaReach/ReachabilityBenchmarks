using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

# ==============================================================================
# Decomposition-based approach for the ISS model
# ==============================================================================
include("iss_BFFPSV18.jl")

import IntervalArithmetic
const IT = Interval{Float64, IntervalArithmetic.Interval{Float64}}

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# tune parameters
tune!(SUITE)

# run the benchmarks
results = run(SUITE, verbose=true)

# return the sample with the smallest time value in each test
println("minimum: ", minimum(results))

# return the median for each test
println("median: ", median(results))

# ==============================================================================
# Create plots
# ==============================================================================

function plot_dashed_lines!(v, time_horizon=20.)
    plot!(x->x, x->v, 0.0, time_horizon, line=2, color="red", linestyle=:dash,
          legend=nothing)
    plot!(x->x, x->-v, 0.0, time_horizon, line=2, color="red", linestyle=:dash,
          legend=nothing)
end

projection_matrix = C[:, 136:270]

𝑂_ISS01[:mode] = "reach"
𝑂_ISS01[:projection_matrix] = projection_matrix
opC = BFFPSV18(merge(𝑂_dense_ISS01, Options(:block_options_iter => nothing)))
res = solve(iss_TV, 𝑂_ISS01, op=opC)

res = ReachSolution(
    [ReachSet(Interval(rs.t_start, rs.t_end) × rs.X, rs.t_start, rs.t_end)
        for rs in res.Xk],
    Options(:plot_vars => [0, 1]))

plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"y_3\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[0., 10., 20.], ytick=[-7e-4, -5e-4, 0., 5e-4, 7e-4],
     xlims=(0., 20.), ylims=(-8e-4, 8e-4),
     title=L"\textcolor{white}{.}\cdot 10^{-4}", title_location=:left,
     titlefont=font(35, "Times"),
     yformatter = yi -> yi == 0 ? "0" : "$(Int(yi*1e4))",
     bottom_margin=6mm, left_margin=2mm, right_margin=2mm,
     size=(1000, 1000))
plot_dashed_lines!(0.0005)
plot_dashed_lines!(0.0007)
savefig(@relpath "ISSF01.png")

𝑂_ISS02[:mode] = "reach"
𝑂_ISS02[:projection_matrix] = projection_matrix
opC = BFFPSV18(merge(𝑂_dense_ISS02, Options(:block_options_iter => nothing)))
res = solve(iss_CONST, 𝑂_ISS02, op=opC)

res = ReachSolution(
    [ReachSet(Interval(rs.t_start, rs.t_end) × rs.X, rs.t_start, rs.t_end)
        for rs in res.Xk],
    Options(:plot_vars => [0, 1]))

plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"y_3\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[0., 10., 20.], ytick=[-17e-5, 0., 17e-5],
     xlims=(0., 20.), ylims=(-2e-4, 2e-4),
     title=L"\textcolor{white}{.}\cdot 10^{-4}", title_location=:left,
     titlefont=font(35, "Times"),
     yformatter = yi -> yi == 0 ? "0" : "$(round(yi*1e4, digits=1))",
     bottom_margin=6mm, left_margin=2mm, right_margin=2mm,
     size=(1000, 1000))
plot_dashed_lines!(0.00017)
savefig(@relpath "ISSC01.png")
