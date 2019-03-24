using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["Quadrotor"] = BenchmarkGroup()

# ==============================================================================
# Jet-based approach using Taylor Models
# ==============================================================================
include("quad.jl")

# benchmark settings
𝑂 = Options(:t0=>0.0, :T=>5.0, :abs_tol=>1e-7, :orderT=>5, :orderQ=>2,
            :max_steps=>500, :property=>quad_property)

# first run
tTM, xTM = quad_TMJets(; t0=𝑂[:t0], T=𝑂[:T], abs_tol=𝑂[:abs_tol],
                orderT=𝑂[:orderT], orderQ=𝑂[:orderQ],
                maxsteps=𝑂[:max_steps], property=𝑂[:property])

# verify that specification holds
@assert all([quad_property(tTM[ind], xTM[ind]) for ind in eachindex(xTM[:])])

# benchmark
SUITE["Quadrotor"]["control property"] = @benchmarkable quad_TMJets(; t0=$𝑂[:t0], T=$𝑂[:T],
                abs_tol=$𝑂[:abs_tol], orderT=$𝑂[:orderT], orderQ=$𝑂[:orderQ],
                maxsteps=$𝑂[:max_steps], property=$𝑂[:property])

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
