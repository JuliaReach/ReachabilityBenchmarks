# =====================================
# MNA 5
#
# system type: affine continuous system
# state dimension: 10913
# =====================================
using ReachabilityBenchmarks, MathematicalSystems, MAT, SparseArrays

function mna5_model()
    file = matopen(@relpath "mna5.mat")

    # system matrix
    A = sparse(read(file, "A"))

    # affine term
    b = sparsevec(19:27, [fill(-0.1, 5); fill(-0.2, 4)], 10913)

    # affine continuous system
    S = @system(x' = Ax + b)

    return S
end
