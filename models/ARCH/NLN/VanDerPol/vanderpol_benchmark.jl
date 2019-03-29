using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["VanDerPol"] = BenchmarkGroup()

# ==============================================================================
# Jet-based approach using Taylor Models
# ==============================================================================
include("vanderpol.jl")

# ----------------------------------------
# Case 1: μ = 1
# ----------------------------------------

# benchmark settings
𝑃, 𝑂 = vanderpol(μ=1)

# algorithm-specific options
𝑂jets = Options(:abs_tol=>1e-10, :orderT=>10, :orderQ=>2, :max_steps=>500)

# first run
sol_1 = solve(𝑃, 𝑂, op=TMJets(𝑂jets))

# verify that specification holds
@assert all([ρ([0.0, 1.0], sol_1.Xk[i].X) < 2.75 for i in eachindex(sol_1.Xk)])

# benchmark
SUITE["VanDerPol"]["μ = 1: x[2] <= 2.75"] = @benchmarkable solve($𝑃, $𝑂, op=TMJets($𝑂jets))

# ----------------------------------------
# Case 2: μ = 2
# ----------------------------------------

# benchmark settings
𝑃, 𝑂 = vanderpol(μ=2, T=8.0,
                 X0=Hyperrectangle(low=[1.55, 2.35], high=[1.85, 2.45]),
                 property=(t,x) -> x[2] < 4.0)

# algorithm-specific options
𝑂jets = Options(:abs_tol=>1e-10, :orderT=>10, :orderQ=>2, :max_steps=>500)

# first run
sol_2 = solve(𝑃, 𝑂, op=TMJets(𝑂jets))

# verify that specification holds
@assert all([ρ([0.0, 1.0], sol_2.Xk[i].X) < 4.0 for i in eachindex(sol_2.Xk)])

# benchmark
SUITE["VanDerPol"]["μ = 2: x[2] <= 4.0"] = @benchmarkable solve($𝑃, $𝑂, op=TMJets($𝑂jets))

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

# ==============================================================================
# Create plots
# ==============================================================================

plot(sol_1,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"x_{1}\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_{2}\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[-3., -2., -1., 0., 1., 2., 3.], ytick=[-3., -2., -1., 0., 1., 2., 3.],
     xlims=(-3., 3.), ylims=(-3., 3.),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="blue")

plot!(x->x, x->2.75, -3., 3., line=2, color="red", linestyle=:dash, legend=nothing)
savefig(@relpath "vanderpol_case_1.png")

plot(sol_2,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"x_{1}\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_{2}\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[-3., -2., -1., 0., 1., 2., 3.], ytick=[-3., -2., -1., 0., 1., 2., 3.],
     xlims=(-3., 3.), ylims=(-3., 3.),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="blue")

plot!(x->x, x->4.0, -3., 3., line=2, color="red", linestyle=:dash, legend=nothing)
savefig(@relpath "vanderpol_case_2.png")
