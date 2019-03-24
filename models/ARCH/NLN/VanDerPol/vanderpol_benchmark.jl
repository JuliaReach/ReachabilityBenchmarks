using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["VanDerPol"] = BenchmarkGroup()

# ==============================================================================
# Jet-based approach using Taylor Models
# ==============================================================================
include("vanderpol.jl")

# benchmark settings
time_horizon = 7.0
𝑂 = Options(:T=>time_horizon, :mode=>"check", :property=>(t, x) -> x[2] < 2.75)

# algorithm-specific options
𝑂jets = Options(:abs_tol=>1e-10, :orderT=>10, :orderQ=>2, :max_steps=>500)

# first run
sol = solve(𝑃, 𝑂, op=TMJets(𝑂jets))

# verify that specification holds
@assert all([ρ([0.0, 1.0], sol.Xk[i].X) < 2.75 for i in eachindex(sol.Xk)])

# benchmark
SUITE["VanDerPol"]["x[2] <= 2.75"] = @benchmarkable solve($𝑃, $𝑂, op=TMJets($𝑂jets))

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
     xlab=L"x_{1}\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_{2}\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[-3., -2., -1., 0., 1., 2., 3.], ytick=[-3., -2., -1., 0., 1., 2., 3.],
     xlims=(-3., 3.), ylims=(-3., 3.),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="lightblue")

plot!(x->x, x->2.75, -3., 3., line=2, color="red", linestyle=:dash, legend=nothing)
savefig(@relpath "vanderpol.png")
