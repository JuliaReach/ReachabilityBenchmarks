# =============================================================================
# To recreate Figure 2, run the following code from the REPL.
#
# julia> include("create_figure_2.jl")
#
# By default, this script requires that you have installed the GR Plots backend.
# See create_figure_5.jl for recommended setups of other plotting backends.
# =============================================================================

using LazySets, LazySets.Approximations, Random
if VERSION >= v"0.7"
    using SparseArrays
else
    using Compat
end

include("plotting.jl")

# set random seed
rng = Random.GLOBAL_RNG
seed = 1234
@static if VERSION < v"0.7-"
    return Random.srand(rng, seed)
else
    return Random.seed!(rng, seed)
end

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

plot(res[1], color=:lightblue, opacity=0.5,
             tickfont=font(20, "Times"), guidefontsize=24,
             xlab=L"x_1\raisebox{-1mm}{\textcolor{white}{.}}",
             ylab=L"x_{50}\raisebox{2mm}{\textcolor{white}{.}}",
             xtick=[0.4, 0.8, 1.2], ytick=[0.8, 1.0, 1.2],
             bottom_margin=6mm, left_margin=6mm)
plot!(res[2], color=:green, opacity=0.5)
plot!(res[3], color=:red, opacity=0.8)

savefig("Figure 2.pdf")
