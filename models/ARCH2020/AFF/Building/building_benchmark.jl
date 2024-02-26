using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["Build"] = BenchmarkGroup()

# ==============================================================================
# Decomposition-based approach
# ==============================================================================
include("building_GLGM06.jl")

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# # tune parameters
# tune!(SUITE)
# 
# # run the benchmarks
# results = run(SUITE, verbose=false)
# 
# # return the sample with the smallest time value in each test
# println("minimum: ", minimum(results))
# 
# # return the median for each test
# println("median: ", median(results))

# ==============================================================================
# Create plots
# ==============================================================================

function plot_dashed_line!(v, time_horizon=20.0)
    return plot!(x -> x, x -> v, 0.0, time_horizon; line=2, color="red", linestyle=:dash,
                 legend=nothing)
end

plot_vars = (0, 25)

res = solve(build_TV; alg=algo_dense, T=time_horizon)

plot(res;
     vars=plot_vars,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t",
     ylab=L"x_{25}",
     xtick=[0.0, 0.5, 1.0], ytick=[-6e-3, -3e-3, 0.0, 3e-3, 6e-3],
     xlims=(0.0, 1.0), ylims=(-9e-3, 6e-3),
     bottom_margin=0mm, left_margin=1mm, right_margin=6mm, top_margin=3mm,
     size=(1000, 1000))
plot_dashed_line!(4e-3)
plot_dashed_line!(5.1e-3)
savefig(@relpath "BLDF01_time_horizon_1.png")

# plot(res,
#      vars=plot_vars,
#      tickfont=font(30, "Times"), guidefontsize=45,
#      xlab=L"t",
#      ylab=L"x_{25}",
#      xtick=[0., 10., 20.], ytick=[-6e-3, -3e-3, 0., 3e-3, 6e-3],
#      bottom_margin=6mm, left_margin=2mm, top_margin=3mm,
#      ylims=(-Inf, 6e-3),
#      size=(1000, 1000))
# plot_dashed_line!(4e-3)
# plot_dashed_line!(5.1e-3)
# savefig(@relpath "BLDF01_time_horizon_20.png")
# 
# res = solve(build_CONST, algo=algo_dense, T=time_horizon)
# 
# plot(res,
#      vars=plot_vars,
#      tickfont=font(30, "Times"), guidefontsize=45,
#      xlab=L"t",
#      ylab=L"x_{25}",
#      xtick=[0., 0.5, 1.], ytick=[-6e-3, -3e-3, 0., 3e-3, 6e-3],
#      xlims=(0., 1.), ylims=(-Inf, 6e-3),
#      bottom_margin=0mm, left_margin=1mm, right_margin=6mm, top_margin=3mm,
#      size=(1000, 1000))
# plot_dashed_line!(4e-3)
# plot_dashed_line!(5.1e-3)
# savefig(@relpath "BLDC01_time_horizon_1.png")
# 
# plot(res,
#      vars=plot_vars,
#      tickfont=font(30, "Times"), guidefontsize=45,
#      xlab=L"t",
#      ylab=L"x_{25}",
#      xtick=[0., 10., 20.], ytick=[-6e-3, -3e-3, 0., 3e-3, 6e-3],
#      bottom_margin=6mm, left_margin=2mm, top_margin=3mm,
#      ylims=(-Inf, 6e-3),
#      size=(1000, 1000))
# plot_dashed_line!(4e-3)
# plot_dashed_line!(5.1e-3)
# savefig(@relpath "BLDC01_time_horizon_20.png")
