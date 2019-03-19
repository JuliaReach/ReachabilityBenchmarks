using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
SUITE["Spacecraft"] = BenchmarkGroup()

include("SpacecraftRendezvous.jl")

SRNA01, options = spacecraft(abort_time=-1.)
SRA01, _ = spacecraft(abort_time=120.)
SRA02, _ = spacecraft(abort_time=[120., 125.])
SRA03, _ = spacecraft(abort_time=[120., 145.])
SRA04, _ = spacecraft(abort_time=240.)
SRA05, _ = spacecraft(abort_time=[235., 240.])
SRA06, _ = spacecraft(abort_time=[230., 240.])
SRA07, _ = spacecraft(abort_time=[50., 150.])
SRA08, _ = spacecraft(abort_time=[0., 240.])
SRU01, _ = spacecraft(abort_time=260.)
SRU02, _ = spacecraft(abort_time=[0., 260.])

# ==============================================================================
# Decomposition-based approach
# ==============================================================================

# algorithm-specific options
options[:mode] = "check"

𝑂_dense = Options(:partition => [1:5], :δ => 0.04)
𝑂_discrete = merge(𝑂_dense, Options(:discretization => "nobloating", :δ => 0.1))
opC_dense = BFFPSV18(𝑂_dense)
opC_discrete = BFFPSV18(𝑂_discrete)
opD = LazyDiscretePost(:lazy_R⋂I => true, :lazy_R⋂G => true)

# single run to verify that specification holds
res = solve(SRNA01, options, opC_dense, opD)
@assert res.satisfied
res = solve(SRNA01, options, opC_discrete, opD)
@assert res.satisfied
res = solve(SRA01, options, opC_dense, opD)
@assert res.satisfied
res = solve(SRA01, options, opC_discrete, opD)
@assert res.satisfied
res = solve(SRA02, options, opC_dense, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA02, options, opC_discrete, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA03, options, opC_dense, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA03, options, opC_discrete, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA04, options, opC_dense, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA04, options, opC_discrete, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA05, options, opC_dense, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA05, options, opC_discrete, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA06, options, opC_dense, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA06, options, opC_discrete, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA07, options, opC_dense, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA07, options, opC_discrete, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA08, options, opC_dense, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRA08, options, opC_discrete, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRU01, options, opC_dense, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRU01, options, opC_discrete, opD)
@assert res.satisfied  # TODO NEW SETTING
res = solve(SRU02, options, opC_dense, opD)
@assert !res.satisfied
res = solve(SRU02, options, opC_discrete, opD)
@assert !res.satisfied

# benchmark
SUITE["Spacecraft"]["SRNA01-SR02", "dense"] =
    @benchmarkable solve($SRNA01, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRNA01-SR02", "discrete"] =
    @benchmarkable solve($SRNA01, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRA01-SR02", "dense"] =
    @benchmarkable solve($SRA01, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRA01-SR02", "discrete"] =
    @benchmarkable solve($SRA01, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRA02-SR02", "dense"] =
    @benchmarkable solve($SRA02, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRA02-SR02", "discrete"] =
    @benchmarkable solve($SRA02, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRA03-SR02", "dense"] =
    @benchmarkable solve($SRA03, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRA03-SR02", "discrete"] =
    @benchmarkable solve($SRA03, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRA04-SR02", "dense"] =
    @benchmarkable solve($SRA04, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRA04-SR02", "discrete"] =
    @benchmarkable solve($SRA04, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRA05-SR02", "dense"] =
    @benchmarkable solve($SRA05, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRA05-SR02", "discrete"] =
    @benchmarkable solve($SRA05, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRA06-SR02", "dense"] =
    @benchmarkable solve($SRA06, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRA06-SR02", "discrete"] =
    @benchmarkable solve($SRA06, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRA07-SR02", "dense"] =
    @benchmarkable solve($SRA07, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRA07-SR02", "discrete"] =
    @benchmarkable solve($SRA07, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRA08-SR02", "dense"] =
    @benchmarkable solve($SRA08, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRA08-SR02", "discrete"] =
    @benchmarkable solve($SRA08, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRU01-SR02", "dense"] =
    @benchmarkable solve($SRU01, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRU01-SR02", "discrete"] =
    @benchmarkable solve($SRU01, $options, $opC_discrete, $opD)
SUITE["Spacecraft"]["SRU02-SR02", "dense"] =
    @benchmarkable solve($SRU02, $options, $opC_dense, $opD)
SUITE["Spacecraft"]["SRU02-SR02", "discrete"] =
    @benchmarkable solve($SRU02, $options, $opC_discrete, $opD)

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# tune parameters
tune!(SUITE)

# run the benchmarks
results = run(SUITE, verbose=false)

# return the sample with the smallest time value in each test
println("minimum: ", minimum(results))

# return the median for each test
println("median: ", median(results))

# ==============================================================================
# Create plots
# ==============================================================================

options[:mode] = "reach"
options[:plot_vars] = [1, 2]
options[:project_reachset] = true
res = solve(SRNA01, options, opC_dense, opD)
plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_1\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[-1000., -500., 0.], ytick=[-400., -300., -200., -100., 0.],
     xlims=(-1000., 0.), ylims=(-450., 50.),
     bottom_margin=6mm, left_margin=6mm,
     size=(1000, 1000))
savefig("SRNA01_SR02.png")

res = solve(SRA01, options, opC_dense, opD)
plot(res,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"x_1\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[-1000., -500., 0.], ytick=[-400., -300., -200., -100., 0.],
     xlims=(-1000., 150.), ylims=(-450., 50.),
     bottom_margin=6mm, left_margin=6mm,
     size=(1000, 1000))
savefig("SRA01_SR02.png")
