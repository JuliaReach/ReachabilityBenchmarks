using ReachabilityAnalysis

# adapt to your own path here
include("/home/mforets/.julia/dev/ReachabilityAnalysis/test/models/hybrid/embrake.jl")

LazySets.deactivate_assertions()

function LazySets.diameter(R::ReachabilityAnalysis.AbstractLazyReachSet; vars::Int)
    overapproximate(project(R, vars=vars), Interval) |> set |> diameter
end

println("In all models Tsample=1e-4")
println("We use δ=1e-8 and static=true options.")

function final_diams(sol)
    println("final diameter of current I = ",  diameter(sol[end][end], vars=(1)))
    println("final diameter of position x = ",  diameter(sol[end][end], vars=(2)))
end

println("==== no pv no jitter (ζ=0.0) ====")
prob = embrake_no_pv(ζ=0.0, Tsample=1e-4);

println("GLGM06, order 1")
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("GLGM06, order 2")
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("GLGM06, order 3")
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=1000);
final_diams(sol)
GC.gc()

println("==== no pv with jitter (ζ=[-1e-8, 1e-7]) ====")
prob = embrake_no_pv(ζ=[-1e-8, 1e-7], Tsample=1e-4);

println("GLGM06, order 1")
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("GLGM06, order 2")
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("GLGM06, order 3")
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=GLGM06(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=1000);
final_diams(sol)
GC.gc()

println("==== pv1 no jitter (ζ=0.0) ====")
prob = embrake_pv_1(ζ=0.0, Tsample=1e-4);

println("ASB07, order 1")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("ASB07, order 2")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("ASB07, order 3")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=1000);
final_diams(sol)
GC.gc()

println("==== pv1 with jitter (ζ=[-1e-8, 1e-7]) ====")
prob = embrake_pv_1(ζ=[-1e-8, 1e-7], Tsample=1e-4);

println("ASB07, order 1")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("ASB07, order 2")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("ASB07, order 3")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=1000);
final_diams(sol)
GC.gc()

println("==== pv2 no jitter (ζ=0.0) ====")
prob = embrake_pv_2(ζ=0.0, Tsample=1e-4);

println("ASB07, order 1")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("ASB07, order 2")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("ASB07, order 3")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=1000);
final_diams(sol)
GC.gc()

println("==== pv2 with jitter (ζ=[-1e-8, 1e-7]) ====")
prob = embrake_pv_2(ζ=[-1e-8, 1e-7], Tsample=1e-4);

println("ASB07, order 1")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("ASB07, order 2")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=2, static=true, dim=4, ngens=8), max_jumps=1000);
final_diams(sol)
GC.gc()

println("")
println("ASB07, order 3")
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=2);
GC.gc()
@time sol = solve(prob, alg=ASB07(δ=1e-8, max_order=3, static=true, dim=4, ngens=12), max_jumps=1000);
final_diams(sol)
GC.gc()
