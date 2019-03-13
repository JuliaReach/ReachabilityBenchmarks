include("iss.jl")

# ==============================
# Setup options
# ==============================

# general options
𝑂_iss = Options(:T=>time_horizon, :mode=>"check", :projection_matrix=>C)

# algorithm-specific options
𝑂_dense = Options(:δ=>5e-3, :vars=>136:270, :assume_sparse=>true)
𝑂_dense_improved_accuracy = Options(:δ=>6e-4, :vars=>136:270, :assume_sparse=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270])
𝑂_discrete = Options(:discretization=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true)
𝑂_discrete_improved_accuracy = Options(:discretization=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270])

# ==============================
# ISU01 and ISS01
# ==============================

# specification options
𝑂_ISU01 = merge(𝑂_iss, Options(:property=>ISU01))
𝑂_ISS01 = merge(𝑂_iss, Options(:property=>ISS01))

# single run
sol_ISU01_dense = solve(iss_TV, 𝑂_ISU01, op=BFFPSV18(𝑂_dense))
sol_ISS01_dense = solve(iss_TV, 𝑂_ISS01, op=BFFPSV18(𝑂_dense_improved_accuracy))
sol_ISU01_discrete = solve(iss_TV, 𝑂_ISU01, op=BFFPSV18(𝑂_discrete))
sol_ISS01_discrete = solve(iss_TV, 𝑂_ISS01, op=BFFPSV18(𝑂_discrete_improved_accuracy))

# verify that specifications hold
@assert !sol_ISU01_dense.satisfied
@assert sol_ISS01_dense.satisfied
@assert !sol_ISU01_discrete.satisfied
@assert sol_ISS01_discrete.satisfied

# benchmark
SUITE["ISS"]["ISU01", "dense"] = @benchmarkable solve($iss_TV, $𝑂_ISU01, op=BFFPSV18($𝑂_dense))
SUITE["ISS"]["ISS01", "dense"] = @benchmarkable solve($iss_TV, $𝑂_ISS01, op=BFFPSV18($𝑂_dense_improved_accuracy))
SUITE["ISS"]["ISU01", "discrete"] = @benchmarkable solve($iss_TV, $𝑂_ISU01, op=BFFPSV18($𝑂_discrete))
SUITE["ISS"]["ISS01", "discrete"] = @benchmarkable solve($iss_TV, $𝑂_ISS01, op=BFFPSV18($𝑂_discrete_improved_accuracy))

# ==============================
# ISU02 and ISS02
# ==============================

# specification options
𝑂_ISU02 = merge(𝑂_iss, Options(:property=>ISU02))
𝑂_ISS02 = merge(𝑂_iss, Options(:property=>ISS02))

# single run
sol_ISU02_dense = solve(iss_CONST, 𝑂_ISU02, op=BFFPSV18(𝑂_dense))
sol_ISS02_dense = solve(iss_CONST, 𝑂_ISS02, op=BFFPSV18(merge(𝑂_dense_improved_accuracy, Options(:δ=>5e-3))))
sol_ISU02_discrete = solve(iss_CONST, 𝑂_ISU02, op=BFFPSV18(𝑂_discrete))
sol_ISS02_discrete = solve(iss_CONST, 𝑂_ISS02, op=BFFPSV18(𝑂_discrete_improved_accuracy))

# verify that specifications hold
@assert !sol_ISU02_dense.satisfied
@assert sol_ISS02_dense.satisfied
@assert !sol_ISU02_discrete.satisfied
@assert sol_ISS02_discrete.satisfied

# benchmark
SUITE["ISS"]["ISU02", "dense"] = @benchmarkable solve($iss_CONST, $𝑂_ISU02, op=BFFPSV18($𝑂_dense))
SUITE["ISS"]["ISS02", "dense"] = @benchmarkable solve($iss_CONST, $𝑂_ISS02, op=BFFPSV18(merge($𝑂_dense_improved_accuracy, Options(:δ=>5e-3))))
SUITE["ISS"]["ISU02", "discrete"] = @benchmarkable solve($iss_CONST, $𝑂_ISU02, op=BFFPSV18($𝑂_discrete))
SUITE["ISS"]["ISS02", "discrete"] = @benchmarkable solve($iss_CONST, $𝑂_ISS02, op=BFFPSV18($𝑂_discrete_improved_accuracy))
