using ReachabilityAnalysis
include("model.jl")

prob = embrake_pv_1(ζ=1e-7, Tsample=1e-4);
solve(prob, alg=ASB07(δ=1e-8, max_order=3), max_jumps=1000, static=false);
