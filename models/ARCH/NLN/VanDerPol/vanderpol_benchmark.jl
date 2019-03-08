using BenchmarkTools
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["VanDerPol"] = BenchmarkGroup()

# ==============================================================================
# Jet-based approach using Taylor Models
# ==============================================================================
include("vanderpol_TMJets.jl")

# benchmark settings
𝑂 = Options(:t0=>0.0, :T=>7.0, :abs_tol=>1e-1, :orderT=>2, :orderQ=>2,
            :maxsteps=>500, :property=>(t, x) -> x[2] < 2.75)

# first run
tTM, xTM = vanderpol_TMJets(; t0=𝑂[:t0], T=𝑂[:T], abs_tol=𝑂[:abs_tol],
                orderT=𝑂[:orderT], orderQ=𝑂[:orderQ],
                maxsteps=𝑂[:maxsteps], property=𝑂[:property])

# verify that specification holds
@assert all([xTM[ind][2] < 2.75 for ind in eachindex(xTM[:])])

# benchmark
SUITE["VanDerPol"]["x[2] <= 2.75"] = @benchmarkable vanderpol_TMJets(; t0=$𝑂[:t0], T=$𝑂[:T],
                abs_tol=$𝑂[:abs_tol], orderT=$𝑂[:orderT], orderQ=$𝑂[:orderQ],
                maxsteps=$𝑂[:maxsteps], property=$𝑂[:property])

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
