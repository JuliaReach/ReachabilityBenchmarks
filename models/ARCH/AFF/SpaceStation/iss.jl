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
problem_TV = InitialValueProblem(S, X0)

# ==============================
# ISU01 and ISS01
# ==============================
ISU01 = LinearConstraintProperty([Clause(LinearConstraint(Cvec, 0.0005)), Clause(LinearConstraint(-Cvec, 0.0005))])
ISS01 = LinearConstraintProperty([Clause(LinearConstraint(Cvec, 0.0007)), Clause(LinearConstraint(-Cvec, 0.0007))])

sol = solve(problem_TV, Options(:T=>time_horizon, :mode=>"check", :property=>ISU01, :projection_matrix=>C), op=BFFPSV18(:δ=>5e-3, :vars=>136:270, :assume_sparse=>true))
SUITE["ISS"]["ISU01", "dense"] = @benchmarkable solve($problem_TV, Options(:T=>$time_horizon, :mode=>"check", :property=>$ISU01, :projection_matrix=>$C), op=BFFPSV18(:δ=>5e-3, :vars=>136:270, :assume_sparse=>true))

sol = solve(problem_TV, Options(:T=>time_horizon, :mode=>"check", :property=>ISS01, :projection_matrix=>C), op=BFFPSV18(:δ=>6e-4, :vars=>136:270, :assume_sparse=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270]))
SUITE["ISS"]["ISS01", "dense"] = @benchmarkable solve(problem_TV, Options(:T=>time_horizon, :mode=>"check", :property=>ISS01, :projection_matrix=>C), op=BFFPSV18(:δ=>6e-4, :vars=>136:270, :assume_sparse=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270]))

sol = solve(problem_TV, Options(:T=>time_horizon, :mode=>"check", :property=>ISU01, :projection_matrix=>C), op=BFFPSV18(:approx_model=>"nobloating", :δ=>5e-3, :vars=> 136:270, :assume_sparse=>true))
SUITE["ISS"]["ISU01", "discrete"] = @benchmarkable solve($problem_TV, Options(:T=>$time_horizon, :mode=>"check", :property=>$ISU01, :projection_matrix=>$C), op=BFFPSV18(:approx_model=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true))

sol = solve(problem_TV, Options(:T=>time_horizon, :mode=>"check", :property=>ISS01, :projection_matrix=>C), op=BFFPSV18(:approx_model=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270]))
SUITE["ISS"]["ISS01", "discrete"] = @benchmarkable solve($problem_TV, Options(:T=>$time_horizon, :mode=>"check", :property=>$ISS01, :projection_matrix=>$C), op=BFFPSV18(:approx_model=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270]))

# ==============================
# Constant input
# ==============================
using Reachability:add_dimension
A = sparse(read(file, "A"))
Aext = add_dimension(A, 3)
Aext[1:270, 271:273] = B
S = LinearContinuousSystem(Aext)
X0 = X0 * U
problem_CONST = InitialValueProblem(S, X0)
C = hcat(C, [0.0 0.0 0.0])
Cvec = C[1, :]

# ==============================
# ISU02 and ISS02
# ==============================
ISU02 = LinearConstraintProperty([Clause(LinearConstraint(Cvec, 0.00017)), Clause(LinearConstraint(-Cvec, 0.00017))])
ISS02 = LinearConstraintProperty([Clause(LinearConstraint(Cvec, 0.0005)), Clause(LinearConstraint(-Cvec, 0.0005))])

sol = solve(problem_CONST, Options(:T=>time_horizon, :mode=>"check", :property=>ISU02, :projection_matrix=>C), op=BFFPSV18(:δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :assume_homogeneous=>true))
SUITE["ISS"]["ISU02", "dense"] = @benchmarkable solve($problem_CONST, Options(:T=>$time_horizon, :mode=>"check", :property=>$ISU02, :projection_matrix=>$C), op=BFFPSV18(:δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :assume_homogeneous=>true))

sol = solve(problem_CONST, Options(:T=>time_horizon, :mode=>"check", :property=>ISS02, :projection_matrix=>C), op=BFFPSV18(:δ=>5e-3, :vars=> 136:270, :assume_sparse=>true, :assume_homogeneous=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270]))
SUITE["ISS"]["ISS02", "dense"] = @benchmarkable solve($problem_CONST, Options(:T=>$time_horizon, :mode=>"check", :property=>$ISS02, :projection_matrix=>$C), op=BFFPSV18(:δ=>5e-3, :vars=> 136:270, :assume_sparse=>true, :assume_homogeneous=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270]))

sol = solve(problem_CONST, Options(:T=>time_horizon, :mode=>"check", :property=>ISU02, :projection_matrix=>C), op=BFFPSV18(:approx_model=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :assume_homogeneous=>true))
SUITE["ISS"]["ISU02", "discrete"] = @benchmarkable solve($problem_CONST, Options(:T=>$time_horizon, :mode=>"check", :property=>$ISU02, :projection_matrix=>$C), op=BFFPSV18(:approx_model=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :assume_homogeneous=>true))

sol = solve(problem_CONST, Options(:T=>time_horizon, :mode=>"check", :property=>ISS02, :projection_matrix=>C), op=BFFPSV18(:approx_model=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :assume_homogeneous=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270]))
SUITE["ISS"]["ISS02", "discrete"] = @benchmarkable solve($problem_CONST, Options(:T=>$time_horizon, :mode=>"check", :property=>$ISS02, :projection_matrix=>$C), op=BFFPSV18(:approx_model=>"nobloating", :δ=>5e-3, :vars=>136:270, :assume_sparse=>true, :assume_homogeneous=>true, :lazy_inputs_interval=>-1, :partition=>[1:135, 136:270]))
