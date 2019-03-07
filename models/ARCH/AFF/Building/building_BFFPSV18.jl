include("building.jl")

# ================
# BLDF01 - BDS01
# ================

# general options
𝑂_BLDF01 = Options(:T=>time_horizon, :mode=>"check", :property => pBDS01)

# algorithm-specific options
𝑂_dense = Options(:δ=>0.004, :vars=>[25], :block_options=>Interval, :lazy_inputs_interval=>-1)
𝑂_discrete = Options(:δ=>0.004, :vars=>[25], :block_options=>Interval, :discretization=>"nobloating", :lazy_inputs_interval=>-1)

# single run
sol_BLDF01_dense = solve(build_TV, 𝑂_BLDF01, op=BFFPSV18(𝑂_dense))
sol_BLDF01_discrete = solve(build_TV, 𝑂_BLDF01, op=BFFPSV18(𝑂_discrete))

# verify that specifications hold
@assert sol_BLDF01_dense.satisfied
@assert sol_BLDF01_discrete.satisfied

# benchmark
SUITE["Build"]["BLDF01-BDS01", "dense"] = @benchmarkable solve($build_TV, $𝑂_BLDF01, op=BFFPSV18($𝑂_dense))
SUITE["Build"]["BLDF01-BDS01", "discrete"] = @benchmarkable solve($build_TV, $𝑂_BLDF01, op=BFFPSV18($𝑂_discrete))

# ================
# BLDC01 - BDS01
# ================

# general options
𝑂_BLDC01 = Options(:T=>time_horizon, :mode=>"check", :property=>pBLDC01)

# algorithm-specific options
𝑂_dense = Options(:δ=>0.004, :vars=>[25], :block_options=>Interval, :lazy_inputs_interval=>-1)
𝑂_discrete = Options(:δ=>0.004, :vars=>[25], :block_options=>Interval, :discretization=>"nobloating", :lazy_inputs_interval=>-1)

# single run
sol_BLDC01_dense = solve(build_CONST, 𝑂_BLDC01, op=BFFPSV18(𝑂_dense))
sol_BLDC01_discrete = solve(build_CONST, 𝑂_BLDC01, op=BFFPSV18(𝑂_discrete))

# verify that specifications hold
@assert sol_BLDC01_dense.satisfied
@assert sol_BLDC01_discrete.satisfied

# register in benchmark suite
SUITE["Build"]["BLDC01-BDS01", "dense"] = @benchmarkable solve($build_CONST, $𝑂_BLDC01, op=BFFPSV18($𝑂_dense))
SUITE["Build"]["BLDC01-BDS01", "discrete"] = @benchmarkable solve($build_CONST, $𝑂_BLDC01, op=BFFPSV18($𝑂_discrete))
