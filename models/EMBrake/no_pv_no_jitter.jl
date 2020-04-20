using ReachabilityAnalysis
include("model.jl")

prob = embrake_no_pv(ζ=0.0, Tsample=1e-4);
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=3, static=false), max_jumps=1000);
