using ReachabilityAnalysis
include("model.jl")

prob = embrake_pv_2(ζ=0.0, Tsample=1e-4, χ=5.0);
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=3, static=false), max_jumps=1000);
