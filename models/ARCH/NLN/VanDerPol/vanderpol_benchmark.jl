using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["VanDerPol"] = BenchmarkGroup()

# ==============================================================================
# Jet-based approach using Taylor Models
# ==============================================================================
include("vanderpol_TMJets.jl")

# benchmark settings
𝑂 = Options(:T=>7.0, :mode=>"check", :property=>(t, x) -> x[2] < 2.75)

# algorithm-specific options
𝑂jets = Options(:abs_tol=>1e-1, :orderT=>2, :orderQ=>2, :maxsteps=>500)

# first run
sol = solve(𝑃, 𝑂, op=TMJets(𝑂jets))

# verify that specification holds
all([ρ([0.0, 1.0], sol.Xk[i].X) < 2.75 for i in eachindex(sol.Xk)])

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
