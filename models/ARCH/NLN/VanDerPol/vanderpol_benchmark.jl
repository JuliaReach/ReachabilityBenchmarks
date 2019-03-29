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

X0_μ2 = Hyperrectangle(low=[1.55, 2.35], high=[1.85, 2.45])

# algorithm-specific options
𝑂jets = Options(:abs_tol=>1e-10, :orderT=>8, :orderQ=>2, :max_steps=>500)

# the idea is to split the initial states along the "x" direction
# n controls the splitting
nsplits_x = 8

function compute_μ2(;n::Int=8, validate=true)
    sol_2 = []
    for X0i in split(X0_μ2, n, 1)
        𝑃, 𝑂 = vanderpol(μ=2.0, T=8.0, X0=X0i, property=(t,x) -> x[2] < 4.0)
        𝑂jets = Options(:abs_tol=>1e-10, :orderT=>8, :orderQ=>1, :max_steps=>500)
        sol_2_i = solve(𝑃, 𝑂, op=TMJets(𝑂jets))
        push!(sol_2, sol_2_i)
        if validate
            # verify that specification holds
            @assert all([ρ([0.0, 1.0], sol_2_i.Xk[i].X) < 4.0 for i in eachindex(sol_2_i.Xk)])
        end
    end
    return sol_2
end

sol_2 = compute_μ2(n=nsplits_x, validate=true)

# benchmark
SUITE["VanDerPol"]["μ = 2: x[2] <= 4.0"] = @benchmarkable compute_μ2(n=$nsplits_x, validate=false)

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

# --------------------------
# Case 1
# --------------------------

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

# --------------------------
# Case 2
# --------------------------

plot_2 = plot(x->x, x->4.0, -2.5, 3., line=2, color="red", linestyle=:dash, legend=nothing)

for i in 1:nsplits_x
    plot!(plot_2, sol_2[i], tickfont=font(30, "Times"), guidefontsize=45,
                   xlab=L"x_{1}\raisebox{-0.5mm}{\textcolor{white}{.}}",
                   ylab=L"x_{2}\raisebox{2mm}{\textcolor{white}{.}}",
                   xtick=[-2., -1., 0., 1., 2., 3.], ytick=[-4., -3., -2., -1., 0., 1., 2., 3., 4.],
                   xlims=(-2.5, 3.), ylims=(-4.5, 4.),
                   bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
                   size=(1000, 1000), color="red", linewidth=0.0, linecolor="red", alpha=.5)
end

savefig(plot_2, @relpath "vanderpol_case_2.png")

# --------------------------
# Cases 1 and 2 overlapped
# --------------------------

plot_all = plot(x->x, x->2.75, -2.5, 3., line=2, color="red", linestyle=:dash, legend=nothing)
plot!(plot_all, x->x, x->4.0, -2.5, 3., line=2, color="red", linestyle=:dash, legend=nothing)

for i in 1:nsplits_x
    plot!(plot_all, sol_2[i], tickfont=font(30, "Times"), guidefontsize=45,
                   xlab=L"x_{1}\raisebox{-0.5mm}{\textcolor{white}{.}}",
                   ylab=L"x_{2}\raisebox{2mm}{\textcolor{white}{.}}",
                   xtick=[-2., -1., 0., 1., 2., 3.], ytick=[-4., -3., -2., -1., 0., 1., 2., 3., 4.],
                   xlims=(-2.5, 3.), ylims=(-4.5, 4.),
                   bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
                   size=(1000, 1000), color="red", linewidth=0.0, linecolor="red", alpha=.5)
end

plot!(sol_1, tickfont=font(30, "Times"), guidefontsize=45,
                       xlab=L"x_{1}\raisebox{-0.5mm}{\textcolor{white}{.}}",
                       ylab=L"x_{2}\raisebox{2mm}{\textcolor{white}{.}}",
                       xtick=[-2., -1., 0., 1., 2., 3.], ytick=[-3., -2., -1., 0., 1., 2., 3., 4.],
                       xlims=(-2.5, 3.), ylims=(-4.5, 4.),
                       bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
                       size=(1000, 1000), color="blue", linewidth=1., linecolor="blue", alpha=.8)

savefig(plot_all, "vanderpol_case_all.png")
