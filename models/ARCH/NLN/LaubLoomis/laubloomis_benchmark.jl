using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["LaubLoomis"] = BenchmarkGroup()

# ==============================================================================
# Jet-based approach using Taylor Models
# ==============================================================================
include("laubloomis.jl")

# --- Case 1: smaller initial states ---
𝑃, 𝑂 = laubloomis(W=0.01, property=(t,x)->x[4] < 4.5)

𝑂₁ = Options(:abs_tol=>1e-10, :orderT=>7, :orderQ=>2, :max_steps=>1000)

# first run
sol = solve(𝑃, 𝑂, op=TMJets(𝑂₁))

# verify that specification holds
v4 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
@assert all([ρ(v4, sol.Xk[i].X) < 4.5 for i in eachindex(sol.Xk)])

# benchmark
SUITE["LaubLoomis"]["W=0.01"] = @benchmarkable solve($𝑃, $𝑂, op=TMJets($𝑂₁))

# --- Case 2: larger initial states ---
𝑃, 𝑂 = laubloomis(W=0.1, property=(t,x)->x[4] < 5.0)

𝑂₂ = copy(𝑂₁)

# first run
sol = solve(𝑃, 𝑂, op=TMJets(𝑂₂))

# verify that specification holds
@assert all([ρ(v4, sol.Xk[i].X) < 5.0 for i in eachindex(sol.Xk)])

# benchmark
SUITE["LaubLoomis"]["W=0.1"] = @benchmarkable solve($𝑃, $𝑂, op=TMJets($𝑂₂))

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
