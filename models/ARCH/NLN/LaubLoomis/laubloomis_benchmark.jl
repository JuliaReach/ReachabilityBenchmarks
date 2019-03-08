using BenchmarkTools
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["LaubLoomis"] = BenchmarkGroup()

# ==============================================================================
# Jet-based approach using Taylor Models
# ==============================================================================
include("laubloomis_TMJets.jl")

# --- Case 1: smaller initial states ---

𝑂₁ = Options(:t0=>0.0, :T=>20.0, :W=>0.01, :abs_tol=>1e-18,
              :orderT=>12, :orderQ=>2, :maxsteps=>1000, :property=>(t,x)->x[4] < 4.5)

# first run
tTM, xTM = laubloomis_TMJets(; t0=𝑂₁[:t0], T=𝑂₁[:T], W=𝑂₁[:W],
                abs_tol=𝑂₁[:abs_tol], orderT=𝑂₁[:orderT], orderQ=𝑂₁[:orderQ],
                maxsteps=𝑂₁[:maxsteps], property=𝑂₁[:property])

# verify that specification holds
@assert all([xTM[ind][4] < 4.5 for ind in eachindex(xTM[:])])

# benchmark
SUITE["LaubLoomis"]["W=0.01"] = @benchmarkable laubloomis_TMJets(; t0=$𝑂₁[:t0], T=$𝑂₁[:T], W=$𝑂₁[:W],
                abs_tol=$𝑂₁[:abs_tol], orderT=$𝑂₁[:orderT], orderQ=$𝑂₁[:orderQ],
                maxsteps=$𝑂₁[:maxsteps], property=$𝑂₁[:property])

# --- Case 2: larger initial states ---

𝑂₂ = Options(:t0=>0.0, :T=>20.0, :W=>0.1, :abs_tol=>1e-24,
              :orderT=>16, :orderQ=>2, :maxsteps=>1000, :property=>(t,x)->x[4] < 5.0)

# first run
tTM, xTM = laubloomis_TMJets(; t0=𝑂₂[:t0], T=𝑂₂[:T], W=𝑂₂[:W],
                abs_tol=𝑂₂[:abs_tol], orderT=𝑂₂[:orderT], orderQ=𝑂₂[:orderQ],
                maxsteps=𝑂₂[:maxsteps], property=𝑂₂[:property])

# verify that specification holds
@assert all([xTM[ind][4] < 5.0 for ind in eachindex(xTM[:])])

# verify tighter specification
#@assert all([xTM[ind][4] < 4.4 for ind in eachindex(xTM[:])])

# benchmark
SUITE["LaubLoomis"]["W=0.1"] = @benchmarkable tTM, xTM = laubloomis_TMJets(; t0=$𝑂₂[:t0], T=$𝑂₂[:T], W=$𝑂₂[:W],
                abs_tol=$𝑂₂[:abs_tol], orderT=$𝑂₂[:orderT], orderQ=$𝑂₂[:orderQ],
                maxsteps=$𝑂₂[:maxsteps], property=$𝑂₂[:property])

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
