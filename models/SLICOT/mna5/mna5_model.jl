# ===========================================================
# MNA 5
#
# system type: LTI system
# state dimension: 10913
# input dimension: 0 (resp. 9 constant inputs)
# ===========================================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT, SparseArrays

function mna5_model()
    file = matopen(@relpath "mna5.mat")

    # system matrix
    A = sparse(read(file, "A"))

    # affine term
    b = sparsevec(19:27, [fill(-0.1, 5); fill(-0.2, 4)], 10913)

    # continuous LTI system
    S = AffineContinuousSystem(A, b)

    return S
end
