# ==============================================================================
# International Space Station model
# See https://easychair.org/publications/open/4cGr
# ==============================================================================

using MAT, Reachability, MathematicalSystems
using SparseArrays, LinearAlgebra, BenchmarkTools

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
S = ConstrainedLinearControlContinuousSystem(
    A, Matrix(1.0I, n, n), nothing, B * U)
iss_TV = InitialValueProblem(S, X0)

# specifications for time-varying input
ISU01 = SafeStatesProperty(HPolyhedron([LinearConstraint(Cvec, 0.0005),
                                        LinearConstraint(-Cvec, 0.0005)]))
ISS01 = SafeStatesProperty(HPolyhedron([LinearConstraint(Cvec, 0.0007),
                                        LinearConstraint(-Cvec, 0.0007)]))

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

# specifications for constant input
ISU02 = SafeStatesProperty(HPolyhedron([LinearConstraint(Cvec, 0.00017),
                                        LinearConstraint(-Cvec, 0.00017)]))
ISS02 = SafeStatesProperty(HPolyhedron([LinearConstraint(Cvec, 0.0005),
                                        LinearConstraint(-Cvec, 0.0005)]))
