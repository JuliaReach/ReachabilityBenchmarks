# ==============================================================================
# Include decomposition-based approach for the ISS model
# ==============================================================================
include("iss_BFFPSV18.jl")

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

#=
Uncomment to re-use parameters file

# From BenchmarkTools.jl
# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `suite` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(SUITE, BenchmarkTools.load(paramspath)[1], :evals);
else
    tune!(SUITE)
    BenchmarkTools.save(paramspath, params(SUITE));
end
=#

results = run(SUITE, verbose=true)

# return the samples with the smallest time value in each test
println(minimum(results))
