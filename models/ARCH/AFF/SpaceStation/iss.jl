using MAT, Reachability, MathematicalSystems
using SparseArrays

SUITE = BenchmarkGroup()

SUITE["ISS"] = BenchmarkGroup()

# ==============================
# Load model
# ==============================
file = matopen(@relpath "iss.mat")
A = sparse(read(file, "A"))
B = read(file, "B")
C = Matrix(read(file, "C")[3, :]')
n = size(A, 1)
Cvec = C[:]
time_horizon = 20.0
X0 = BallInf(zeros(n), 0.0001)

# ==============================
# Time-varying input
# ==============================
U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.])
S = ConstrainedLinearControlContinuousSystem(A, Matrix(1.0I, n, n), nothing, ConstantInput(B * U))
iss_TV = InitialValueProblem(S, X0)

# ==============================
# Setup options
# ==============================

# general options
𝑂_iss = Options(:T=>time_horizon, :mode=>"check", :projection_matrix=>C)

# algorithm-specific options
𝑂_dense = Options(:δ=>5e-3, :vars=>136:270, :assume_sparse=>true)
𝑂_dense_improved_accuracy = Options(:δ=>6e-4, :vars=>136:270, :assume_sparse=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270])
𝑂_discrete = Options(:approx_model=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true)
𝑂_discrete_improved_accuracy = Options(:approx_model=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270])

# ==============================
# ISU01 and ISS01
# ==============================
ISU01 = LinearConstraintProperty([Clause(LinearConstraint(Cvec, 0.0005)), Clause(LinearConstraint(-Cvec, 0.0005))])
ISS01 = LinearConstraintProperty([Clause(LinearConstraint(Cvec, 0.0007)), Clause(LinearConstraint(-Cvec, 0.0007))])

# specification options
𝑂_ISU01 = merge(𝑂_iss, Options(:property=>ISU01))
𝑂_ISS01 = merge(𝑂_iss, Options(:property=>ISS01))

# get things compiled
sol_ISU01_dense = solve(iss_TV, 𝑂_ISU01, op=BFFPSV18(𝑂_dense))
sol_ISS01_dense = solve(iss_TV, 𝑂_ISS01, op=BFFPSV18(𝑂_dense_improved_accuracy))
sol_ISU01_discrete = solve(iss_TV, 𝑂_ISU01, op=BFFPSV18(𝑂_discrete))
sol_ISS01_discrete = solve(iss_TV, 𝑂_ISS01, op=BFFPSV18(𝑂_discrete_improved_accuracy))

# verify that specifications hold
#@assert sol_ISU01_dense.satisfied
#@assert sol_ISS01_dense.satisfied
#@assert sol_ISU01_discrete.satisfied
#@assert sol_ISS01_discrete.satisfied

# benchmark
SUITE["ISS"]["ISU01", "dense"] = @benchmarkable solve($iss_TV, $𝑂_ISU01, op=BFFPSV18($𝑂_dense))
SUITE["ISS"]["ISS01", "dense"] = @benchmarkable solve($iss_TV, $𝑂_ISS01, op=BFFPSV18($𝑂_dense_improved_accuracy))
SUITE["ISS"]["ISU01", "discrete"] = @benchmarkable solve($iss_TV, $𝑂_ISU01, op=BFFPSV18($𝑂_discrete))
SUITE["ISS"]["ISS01", "discrete"] = @benchmarkable solve($iss_TV, $𝑂_ISS01, op=BFFPSV18($𝑂_discrete_improved_accuracy))

# ==============================
# Constant input
# ==============================
using Reachability:add_dimension
A = sparse(read(file, "A"))
Aext = add_dimension(A, 3)
Aext[1:270, 271:273] = B
S = LinearContinuousSystem(Aext)
X0 = X0 * U
iss_CONST = InitialValueProblem(S, X0)
C = hcat(C, [0.0 0.0 0.0])
Cvec = C[1, :]

# ==============================
# ISU02 and ISS02
# ==============================
ISU02 = LinearConstraintProperty([Clause(LinearConstraint(Cvec, 0.00017)), Clause(LinearConstraint(-Cvec, 0.00017))])
ISS02 = LinearConstraintProperty([Clause(LinearConstraint(Cvec, 0.0005)), Clause(LinearConstraint(-Cvec, 0.0005))])

# specification options
𝑂_ISU02 = merge(𝑂_iss, Options(:property=>ISU02))
𝑂_ISS02 = merge(𝑂_iss, Options(:property=>ISS02))

# get things compiled
sol_ISU02_dense = solve(iss_CONST, 𝑂_ISU02, op=BFFPSV18(𝑂_dense))
sol_ISS02_dense = solve(iss_CONST, 𝑂_ISS02, op=BFFPSV18(merge(𝑂_dense_improved_accuracy, Options(:δ=>5e-3))))
sol_ISU02_discrete = solve(iss_CONST, 𝑂_ISU02, op=BFFPSV18(𝑂_discrete))
sol_ISS02_discrete = solve(iss_CONST, 𝑂_ISS02, op=BFFPSV18(𝑂_discrete_improved_accuracy))

# verify that specifications hold
#@assert sol_ISU02_dense.satisfied
#@assert sol_ISS02_dense.satisfied
#@assert sol_ISU02_discrete.satisfied
#@assert sol_ISS02_discrete.satisfied

# benchmark
SUITE["ISS"]["ISU02", "dense"] = @benchmarkable solve($iss_CONST, $𝑂_ISU02, op=BFFPSV18($𝑂_dense))
SUITE["ISS"]["ISS02", "dense"] = @benchmarkable solve($iss_CONST, $𝑂_ISS02, op=BFFPSV18(merge($𝑂_dense_improved_accuracy, Options(:δ=>5e-3))))
SUITE["ISS"]["ISU02", "discrete"] = @benchmarkable solve($iss_CONST, $𝑂_ISU02, op=BFFPSV18($𝑂_discrete))
SUITE["ISS"]["ISS02", "discrete"] = @benchmarkable solve($iss_CONST, $𝑂_ISS02, op=BFFPSV18($𝑂_discrete_improved_accuracy))
