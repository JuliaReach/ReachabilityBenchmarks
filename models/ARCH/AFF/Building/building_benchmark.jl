using BenchmarkTools
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["Build"] = BenchmarkGroup()

# ==============================================================================
# Decomposition-based approach for the Building model
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
