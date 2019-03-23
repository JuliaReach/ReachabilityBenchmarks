using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["Build"] = BenchmarkGroup()

# ==============================================================================
# Decomposition-based approach
# ==============================================================================
include("building_BFFPSV18.jl")

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# tune parameters
tune!(SUITE)

# run the benchmarks
results = run(SUITE, verbose=false)

# return the sample with the smallest time value in each test
println("minimum: ", minimum(results))

# return the median for each test
println("median: ", median(results))

# ==============================================================================
# Create plots
# ==============================================================================

function plot_dashed_line!(v, time_horizon=20.)
    plot!(x->x, x->v, 0.0, time_horizon, line=2, color="red", linestyle=:dash,
          legend=nothing)
end

𝑂_BLDF01[:mode] = "reach"
𝑂_BLDF01[:plot_vars] = [0, 25]
𝑂_BLDF01[:project_reachset] = true
res = solve(build_TV, 𝑂_BLDF01, op=BFFPSV18(𝑂_dense_BLDF01))

plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_{25}\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[0., 0.5, 1.], ytick=[-6e-3, -3e-3, 0., 3e-3, 6e-3],
     xlims=(0., 1.), ylims=(-Inf, 6e-3),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000))
plot_dashed_line!(4e-3)
plot_dashed_line!(5.1e-3)
savefig("BLDF01_time_horizon_1.png")

plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_{25}\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[0., 10., 20.], ytick=[-6e-3, -3e-3, 0., 3e-3, 6e-3],
     bottom_margin=6mm, left_margin=2mm, top_margin=3mm,
     ylims=(-Inf, 6e-3),
     size=(1000, 1000))
plot_dashed_line!(4e-3)
plot_dashed_line!(5.1e-3)
savefig("BLDF01_time_horizon_20.png")

𝑂_BLDC01[:mode] = "reach"
𝑂_BLDC01[:plot_vars] = [0, 25]
𝑂_BLDC01[:project_reachset] = true
res = solve(build_CONST, 𝑂_BLDC01, op=BFFPSV18(𝑂_dense_BLDC01))

plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_{25}\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[0., 0.5, 1.], ytick=[-6e-3, -3e-3, 0., 3e-3, 6e-3],
     xlims=(0., 1.), ylims=(-Inf, 6e-3),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000))
plot_dashed_line!(4e-3)
plot_dashed_line!(5.1e-3)
savefig("BLDC01_time_horizon_1.png")

plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_{25}\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[0., 10., 20.], ytick=[-6e-3, -3e-3, 0., 3e-3, 6e-3],
     bottom_margin=6mm, left_margin=2mm, top_margin=3mm,
     ylims=(-Inf, 6e-3),
     size=(1000, 1000))
plot_dashed_line!(4e-3)
plot_dashed_line!(5.1e-3)
savefig("BLDC01_time_horizon_20.png")
