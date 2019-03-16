# ==============================================================================
# Building model
# See https://easychair.org/publications/open/4cGr
# ==============================================================================

using MAT, Reachability, MathematicalSystems, SX
using SparseArrays, LinearAlgebra

# ==============================
# Load model
# ==============================
file = @relpath "SX/Building_more_decimals.xml"
H = readsxmodel(file, ST=ConstrainedLinearControlContinuousSystem)

# ===================
# Time-varying input
# ===================

n = size(H.modes[1].A, 1)-1 # the sx model has "time" as a state variable
@assert n == 48
A = H.modes[1].A[1:n, 1:n] 
B = H.modes[1].B[1:n, 1]
U = Hyperrectangle(low=[0.8], high=[1.0])
X = nothing
S = ConstrainedLinearControlContinuousSystem(A, B, X, U)

# specify initial states
center_X0 = [fill(0.000225, 10); fill(0.0, 38)]
radius_X0 = [fill(0.000025, 10); fill(0.0, 14); 0.0001; fill(0.0, 23)]
X0 = Hyperrectangle(center_X0, radius_X0)
build_TV = InitialValueProblem(S, X0)

# specifications
pBDS01 = SafeStatesProperty(
    LinearConstraint(sparsevec([25], [1.0], n), 0.0051)  # x25 <= 0.0051
    )

time_horizon = 20.0

# ===================
# Constant input
# ===================

A = Reachability.add_dimension(A) # add an extra zero row and column
S = LinearContinuousSystem(A)
X0 = X0 × U
build_CONST = InitialValueProblem(S, X0)

# specifications
pBLDC01 = SafeStatesProperty(
    LinearConstraint(sparsevec([25], [1.0], n+1), 0.0051)  # x25 <= 0.0051
    )
