using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["Quadrotor"] = BenchmarkGroup()

# ==============================================================================
# Jet-based approach using Taylor Models
# ==============================================================================
include("quadrotor.jl")

𝑃, 𝑂 = quad()

# algorithm-specific options
𝑂jets = Options(:abs_tol=>1e-7, :orderT=>5, :orderQ=>2, :max_steps=>500)

# first run
sol = solve(𝑃, 𝑂, op=TMJets(𝑂jets))

# benchmark
SUITE["Quadrotor"]["Specification"] = @benchmarkable solve($𝑃, $𝑂, op=TMJets($𝑂jets))

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# tune parameters
tune!(SUITE)

# run the benchmarks
results = run(SUITE, verbose=true)

# return the sample with the smallest time value in each test
println("minimum time for each benchmark:\n", minimum(results))

# return the median for each test
println("median time for each benchmark:\n", median(results))

# ==============================================================================
# Create plots
# ==============================================================================

plot(sol,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_{3}\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[0., 1, 2, 3, 4, 5], ytick=[0.5, 0., 0.5, 1.0, 1.5],
     xlims=(0., 5.), ylims=(-0.5, 1.5),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="blue")

plot!(x->x, x->1.4, 0., 5., line=2, color="red", linestyle=:dash, legend=nothing)
plot!(x->x, x->0.98, 0., 5., line=2, color="red", linestyle=:dash, legend=nothing)
plot!(x->x, x->1.02, 0., 5., line=2, color="red", linestyle=:dash, legend=nothing)
savefig(@relpath "quadrotor.png")
