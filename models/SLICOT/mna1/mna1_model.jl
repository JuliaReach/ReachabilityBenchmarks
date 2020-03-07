# =====================================
# MNA 1
#
# system type: affine continuous system
# state dimension: 578
# =====================================
using ReachabilityBenchmarks, MathematicalSystems, MAT, SparseArrays

function mna1_model()
    file = matopen(@relpath "mna1.mat")

    # system matrix
    A = sparse(read(file, "A"))

    # affine term
    b = sparsevec(570:578, [fill(-0.1, 5); fill(-0.2, 4)], 578)

    # affine continuous system
    S = @system(x' = Ax + b)

    return S
end
