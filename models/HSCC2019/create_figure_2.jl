# To recreate Figure 2, run the following code from the REPL.
# Just type:
#
#     include("create_figure_2.jl")
#
# Note that the code uses random numbers, so the resulting plots and the runtime will likely differ.

using LazySets, LazySets.Approximations, Plots
if VERSION >= v"0.7"
    using SparseArrays
else
    using Compat
end

ENV["GKSwstype"] = "100"

n = 1000; m = 2; δ = 0.1;
A, B = sprandn(n, n, 0.01), randn(n, m);
X0 = BallInf(ones(n), 0.1);
U = Ball2(zeros(m), 1.2);
Y = ConvexHull(SparseMatrixExp(A * δ) * X0 ⊕ δ * B *U, X0);
π = spzeros(2, n) ; π[1, 1] = π[2, 50] = 1;
res = Vector{HPolygon{Float64}}(undef, 3)
εs = [Inf, 0.1, 0.01]
println("warm-up runtimes:")
for (i, ε) in enumerate(εs)
    @time res[i] = overapproximate(π * Y, ε);
end
println("runtimes to obtain an overapproximation from the projection:")
for (i, ε) in enumerate(εs)
    @time res[i] = overapproximate(π * Y, ε);
end

plot(res[1], color=:lightblue, opacity=0.5)
plot!(res[2], color=:green, opacity=0.5)
plot!(res[3], color=:red, opacity=0.5)

# optionally store the picture in a file
savefig("Figure 2.png")
