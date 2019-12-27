# ===========================================================
# MNA 1
#
# system type: LTI system
# state dimension: 578
# input dimension: 0 (resp. 9 constant inputs)
# ===========================================================
using ReachabilityBenchmarks, MathematicalSystems, LazySets, MAT, SparseArrays

function mna1_model()
    file = matopen(@relpath "mna1.mat")

    # system matrix
    A = sparse(read(file, "A"))

    # affine term
    b = sparsevec(570:578, [fill(-0.1, 5); fill(-0.2, 4)], 578)

    # continuous LTI system
    S = AffineContinuousSystem(A, b)

    return S
end
