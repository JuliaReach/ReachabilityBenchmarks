# =====================================
# Projectile
#
# system type: affine continuous system
# state dimension: 4
# =====================================
using MathematicalSystems, SparseArrays

function projectile_model()
    # system matrix
    A = sparse([1, 3], [2, 4], [0.5, 0.7], 4, 4)

    # affine term
    b = sparsevec([4], [-9.81], 4)

    # affine continuous system
    S = @system(x' = Ax + b)

    return S
end
