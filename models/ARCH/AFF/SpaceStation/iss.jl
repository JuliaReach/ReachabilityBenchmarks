# ==============================================================================
# International Space Station model
# See https://easychair.org/publications/open/4cGr
# ==============================================================================

using ReachabilityBenchmarks, MAT, Reachability, MathematicalSystems,
      MathematicalPredicates, SparseArrays

# ==============================
# Load model
# ==============================
file = matopen(@current_path "iss.mat")
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
U = Hyperrectangle(; low=[0.0, 0.8, 0.9], high=[0.1, 1.0, 1.0])
S = ConstrainedLinearControlContinuousSystem(A, B, nothing, U)
iss_TV = InitialValueProblem(S, X0)

# specifications for time-varying input
ISU01 = is_contained_in(HPolyhedron([HalfSpace(Cvec, 0.0005),
                                     HalfSpace(-Cvec, 0.0005)]))
ISS01 = is_contained_in(HPolyhedron([HalfSpace(Cvec, 0.0007),
                                     HalfSpace(-Cvec, 0.0007)]))

# ==============================
# Constant input
# ==============================
using Reachability: add_dimension
A = sparse(read(file, "A"))
Aext = add_dimension(A, 3)
Aext[1:270, 271:273] = B
S = LinearContinuousSystem(Aext)
X0 = X0 * U
iss_CONST = InitialValueProblem(S, X0)
C = hcat(C, [0.0 0.0 0.0])
Cvec = C[1, :]

# specifications for constant input
ISU02 = is_contained_in(HPolyhedron([HalfSpace(Cvec, 0.00017),
                                     HalfSpace(-Cvec, 0.00017)]))
ISS02 = is_contained_in(HPolyhedron([HalfSpace(Cvec, 0.0005),
                                     HalfSpace(-Cvec, 0.0005)]))
