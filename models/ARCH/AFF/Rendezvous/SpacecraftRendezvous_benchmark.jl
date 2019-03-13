using BenchmarkTools
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["Spacecraft"] = BenchmarkGroup()

include("SpacecraftRendezvous.jl")

SR01, options_SR01 = spacecraft(false)
SR02, options_SR02 = spacecraft(true)

# ==============================================================================
# Decomposition-based approach
# ==============================================================================

# algorithm-specific options
options_SR01[:mode] = "check"
options_SR02[:mode] = "check"

𝑂_dense = Options(:partition => [1:5], :δ => 0.04)
𝑂_discrete = merge(𝑂_dense, Options(:discretization => "nobloating", :δ => 0.1))
opC_dense = BFFPSV18(𝑂_dense)
opC_discrete = BFFPSV18(𝑂_discrete)
opD = LazyDiscretePost(:lazy_R⋂I => true, :lazy_R⋂G => true)

# single run to verify that specification holds
res = solve(SR01, options_SR01, opC_dense, opD)
@assert res.satisfied
res = solve(SR01, options_SR01, opC_discrete, opD)
@assert res.satisfied
res = solve(SR02, options_SR02, opC_dense, opD)
@assert res.satisfied
res = solve(SR02, options_SR02, opC_discrete, opD)
@assert res.satisfied

# benchmark
SUITE["Spacecraft"]["SRNA01-SR01", "dense"] =
    @benchmarkable solve($SR01, $options_SR01, $opC_dense, $opD)
SUITE["Spacecraft"]["SRNA01-SR01", "discrete"] =
    @benchmarkable solve($SR01, $options_SR01, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRA01-SR02", "dense"] =
    @benchmarkable solve($SR02, $options_SR02, $opC_dense, $opD)
SUITE["Spacecraft"]["SRA01-SR02", "discrete"] =
    @benchmarkable solve($SR02, $options_SR02, $opC_discrete, $opD)

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
