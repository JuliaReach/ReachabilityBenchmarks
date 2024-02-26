using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["LaubLoomis"] = BenchmarkGroup()

# ==============================================================================
# Jet-based approach using Taylor Models
# ==============================================================================
include("laubloomis.jl")

# ----------------------------------------
# --- Case 1: smaller initial states ---
# ----------------------------------------
𝑃, 𝑂 = laubloomis(; W=0.01, property=(t, x) -> x[4] < 4.5)

𝑂₁ = Options(:abs_tol => 1e-10, :orderT => 7, :orderQ => 1, :max_steps => 1000)

# first run
sol_case_1 = solve(𝑃, 𝑂; op=TMJets(𝑂₁))

# verify that specification holds
v4 = [0.0, 1.0] # the flowpipe has been projected so we check for the second component which is x4
@assert all([ρ(v4, sol_case_1.Xk[i].X) < 4.5 for i in eachindex(sol_case_1.Xk)])

# print width of final box
H = sol_case_1.Xk[end].X
println("width of final box, case W = 0.01 : $(high(H)[2] - low(H)[2])")

# benchmark
SUITE["LaubLoomis"]["W=0.01"] = @benchmarkable solve($𝑃, $𝑂, op=TMJets($𝑂₁))

# ----------------------------------------
# --- Case 2: intermediate initial states ---
# ----------------------------------------
𝑃, 𝑂 = laubloomis(; W=0.05, property=(t, x) -> x[4] < 4.5)

𝑂₂ = Options(:abs_tol => 1e-10, :orderT => 7, :orderQ => 1, :max_steps => 1000)

# first run
sol_case_2 = solve(𝑃, 𝑂; op=TMJets(𝑂₂))

# verify that specification holds
@assert all([ρ(v4, sol_case_2.Xk[i].X) < 4.5 for i in eachindex(sol_case_2.Xk)])

# print width of final box
H = sol_case_2.Xk[end].X
println("width of final box, case W = 0.05 : $(high(H)[2] - low(H)[2])")

# benchmark
SUITE["LaubLoomis"]["W=0.05"] = @benchmarkable solve($𝑃, $𝑂, op=TMJets($𝑂₂))

# ----------------------------------------
# --- Case 3: larger initial states ---
# ----------------------------------------
𝑃, 𝑂 = laubloomis(; W=0.1, property=(t, x) -> x[4] < 5.0)

𝑂₃ = Options(:abs_tol => 1e-10, :orderT => 7, :orderQ => 1, :max_steps => 1000)

# first run
sol_case_3 = solve(𝑃, 𝑂; op=TMJets(𝑂₃))

# verify that specification holds
@assert all([ρ(v4, sol_case_3.Xk[i].X) < 5.0 for i in eachindex(sol_case_3.Xk)])

# print width of final box
H = sol_case_3.Xk[end].X
println("width of final box, case W = 0.1 : $(high(H)[2] - low(H)[2])")

# benchmark
SUITE["LaubLoomis"]["W=0.1"] = @benchmarkable solve($𝑃, $𝑂, op=TMJets($𝑂₃))

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# tune parameters
tune!(SUITE)

# run the benchmarks
results = run(SUITE; verbose=true)

# return the sample with the smallest time value in each test
println("minimum time for each benchmark:\n", minimum(results))

# return the median for each test
println("median time for each benchmark:\n", median(results))

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

plot(sol_case_1;
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{2.0mm}{\textcolor{white}{.}}",
     ylab=L"x_{4}\raisebox{1.2mm}{\textcolor{white}{.}}",
     xtick=[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
     ytick=[2, 2.5, 3, 3.5, 4, 4.5],
     xlims=(0.0, 20.0), ylims=(1.5, 4.5),
     bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="blue")

plot!(x -> x, x -> 4.5, 0.0, 20.0; line=2, color="red", linestyle=:dash, legend=nothing)
savefig(@relpath "laubloomis_case_1.png")

plot(sol_case_2;
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{2.0mm}{\textcolor{white}{.}}",
     ylab=L"x_{4}\raisebox{1.2mm}{\textcolor{white}{.}}",
     xtick=[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
     ytick=[2, 2.5, 3, 3.5, 4, 4.5, 5.0],
     xlims=(0.0, 20.0), ylims=(1.5, 5.0),
     bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="blue")

plot!(x -> x, x -> 5.0, 0.0, 20.0; line=2, color="red", linestyle=:dash, legend=nothing)
savefig(@relpath "laubloomis_case_2.png")

plot(sol_case_3;
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{2.0mm}{\textcolor{white}{.}}",
     ylab=L"x_{4}\raisebox{1.2mm}{\textcolor{white}{.}}",
     xtick=[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
     ytick=[2, 2.5, 3, 3.5, 4, 4.5, 5.0],
     xlims=(0.0, 20.0), ylims=(1.5, 5.0),
     bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="blue")

plot!(x -> x, x -> 5.0, 0.0, 20.0; line=2, color="red", linestyle=:dash, legend=nothing)
savefig(@relpath "laubloomis_case_3.png")

plot(sol_case_1; color="red")

plot!(sol_case_2; alpha=0.6, color="green")

plot!(sol_case_3; alpha=0.2,
      tickfont=font(30, "Times"), guidefontsize=45,
      xlab=L"t\raisebox{2.0mm}{\textcolor{white}{.}}",
      ylab=L"x_{4}\raisebox{1.2mm}{\textcolor{white}{.}}",
      xtick=[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
      ytick=[2, 2.5, 3, 3.5, 4, 4.5, 5.0],
      xlims=(0.0, 20.0), ylims=(1.5, 5.0),
      bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm,
      size=(1000, 1000), color="blue")

plot!(x -> x, x -> 5.0, 0.0, 20.0; line=2, color="red", linestyle=:dash, legend=nothing)
plot!(x -> x, x -> 4.5, 0.0, 20.0; line=2, color="red", linestyle=:dash, legend=nothing)
savefig(@relpath "laubloomis_case_all.png")
