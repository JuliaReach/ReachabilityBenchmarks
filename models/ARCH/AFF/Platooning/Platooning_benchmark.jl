using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["Platooning"] = BenchmarkGroup()

include("Platooning.jl")

PLAD01_BND42, options_PLAD01_BND42 =
    platooning(; deterministic_switching=true, time_horizon=20.,
                 allowed_distance=42.)
PLAD01_BND30, options_PLAD01_BND30 =
    platooning(; deterministic_switching=true, time_horizon=20.,
                 allowed_distance=30.)
PLAN01_UNB50, options_PLAN01_UNB50 =
    platooning(; deterministic_switching=true, time_horizon=Inf,
                 allowed_distance=50.)

# ==============================================================================
# Decomposition-based approach
# ==============================================================================

# algorithm-specific options
options_PLAD01_BND42[:mode] = "check"
options_PLAD01_BND30[:mode] = "check"
options_PLAN01_UNB50[:mode] = "check"

𝑂_common = Options(:partition => [1:10])
𝑂_dense_options_PLAD01_BND42 = merge(𝑂_common, Options(:δ => 0.01))
𝑂_dense_options_PLAD01_BND30 = merge(𝑂_common, Options(:δ => 0.00001))
𝑂_dense_options_PLAN01_UNB50 = merge(𝑂_common, Options(:δ => 0.03))
𝑂_discrete = merge(𝑂_common, Options(:discretization => "nobloating", :δ => 0.05))

opC_dense_PLAD01_BND42 = BFFPSV18(𝑂_dense_options_PLAD01_BND42)
opC_dense_PLAD01_BND30 = BFFPSV18(𝑂_dense_options_PLAD01_BND30)
opC_dense_PLAN01_UNB50 = BFFPSV18(𝑂_dense_options_PLAN01_UNB50)
opC_discrete = BFFPSV18(𝑂_discrete)
opD = LazyDiscretePost(:lazy_R⋂I => true, :lazy_R⋂G => false)

# single run to verify that specification holds
res = solve(PLAD01_BND42, options_PLAD01_BND42, opC_dense_PLAD01_BND42, opD)
@assert res.satisfied
res = solve(PLAD01_BND42, options_PLAD01_BND42, opC_discrete, opD)
@assert res.satisfied
# res = solve(PLAD01_BND30, options_PLAD01_BND30, opC_dense_PLAD01_BND30, opD)
# @assert res.satisfied
# res = solve(PLAD01_BND30, options_PLAD01_BND30, opC_discrete, opD)
# @assert res.satisfied
res = solve(PLAN01_UNB50, options_PLAN01_UNB50, opC_dense_PLAN01_UNB50, opD)
@assert res.satisfied
res = solve(PLAN01_UNB50, options_PLAN01_UNB50, opC_discrete, opD)
@assert res.satisfied

# benchmark
SUITE["Platooning"]["PLAD01_BND42", "dense"] =
    @benchmarkable solve($PLAD01_BND42, $options_PLAD01_BND42, $opC_dense_PLAD01_BND42, $opD)
SUITE["Platooning"]["PLAD01_BND42", "discrete"] =
    @benchmarkable solve($PLAD01_BND42, $options_PLAD01_BND42, $opC_discrete, $opD)
# SUITE["Platooning"]["PLAD01_BND30", "dense"] =
#     @benchmarkable solve($PLAD01_BND30, $options_PLAD01_BND30, $opC_dense_PLAD01_BND30, $opD)
# SUITE["Platooning"]["PLAD01_BND30", "discrete"] =
#     @benchmarkable solve($PLAD01_BND30, $options_PLAD01_BND30, $opC_discrete, $opD)
SUITE["Platooning"]["PLAN01_UNB50", "dense"] =
    @benchmarkable solve($PLAN01_UNB50, $options_PLAN01_UNB50, $opC_dense_PLAN01_UNB50, $opD)
SUITE["Platooning"]["PLAN01_UNB50", "discrete"] =
    @benchmarkable solve($PLAN01_UNB50, $options_PLAN01_UNB50, $opC_discrete, $opD)

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

options_PLAD01_BND42[:mode] = "reach"
options_PLAD01_BND42[:plot_vars] = [0, 1]
options_PLAD01_BND42[:project_reachset] = true
res = solve(PLAD01_BND42, options_PLAD01_BND42, opC_dense_PLAD01_BND42, opD)
plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_1\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[0., 10., 20.], ytick=[-40., -20., 0., 20.],
     bottom_margin=6mm, left_margin=2mm,
     size=(1000, 1000))
plot_dashed_line!(-42.)
savefig("PLAD01_BND42.png")

options_PLAN01_UNB50[:mode] = "reach"
options_PLAN01_UNB50[:plot_vars] = [0, 1]
options_PLAN01_UNB50[:project_reachset] = true
res = solve(PLAN01_UNB50, options_PLAN01_UNB50, opC_dense_PLAN01_UNB50, opD)
time_horizon = res.Xk[end].t_start
plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_1\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[0., 100., 200.], ytick=[-40., -20., 0., 20.],
     bottom_margin=6mm, left_margin=2mm,
     size=(1000, 1000))
plot_dashed_line!(-50., time_horizon)
savefig("PLAN01_UNB50.png")
