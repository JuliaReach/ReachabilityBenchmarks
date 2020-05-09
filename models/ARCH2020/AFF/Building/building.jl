# ==============================================================================
# Building model
# See https://easychair.org/publications/open/4cGr
# ==============================================================================

using ReachabilityBenchmarks, ReachabilityAnalysis, MathematicalSystems, SX,
      MathematicalPredicates, SparseArrays

# ==============================
# Load model
# ==============================
file = @relpath "SpaceEx/Building_more_decimals.xml"
H = readsxmodel(file, ST=ConstrainedLinearControlContinuousSystem)

# ===================
# Time-varying input
# ===================

n = size(H.modes[1].A, 1)-1 # the sx model has "time" as a state variable
@assert n == 48
A = H.modes[1].A[1:n, 1:n] 
B = hcat(H.modes[1].B[1:n, 1])
X = nothing
U = Hyperrectangle(low=[0.8], high=[1.0])
S = ConstrainedLinearControlContinuousSystem(A, B, X, U)

# specify initial states
center_X0 = [fill(0.000225, 10); fill(0.0, 38)]
radius_X0 = [fill(0.000025, 10); fill(0.0, 14); 0.0001; fill(0.0, 23)]
X0 = Hyperrectangle(center_X0, radius_X0)
build_TV = InitialValueProblem(S, X0)

# specifications
pBDS01 = is_contained_in(HalfSpace(sparsevec([25], [1.0], n), 0.0051))  # x25 <= 0.0051

time_horizon = 20.0

# ===================
# Constant input
# ===================

A = ReachabilityAnalysis.add_dimension(A) # add an extra zero row and column
S = LinearContinuousSystem(A)
X0 = X0 × U
build_CONST = InitialValueProblem(S, X0)

# specifications
pBLDC01 = is_contained_in(HalfSpace(sparsevec([25], [1.0], n+1), 0.0051))  # x25 <= 0.0051
