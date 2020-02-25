# ===========================================================================
# Transmission line
#
# system type: LTI system with uncertain dynamics
# state dimension: parametric (2 * η)
# input dimension: 1
#
# Sources:
# [1] M. Althoff, B. H. Krogh, and O. Stursberg. Analyzing Reachability of
#     Linear Dynamic Systems with Parametric Uncertainties. Modeling, Design,
#     and Simulation of Systems with Uncertainties. 2011.
# ===========================================================================
using MathematicalSystems, LazySets, IntervalMatrices, LinearAlgebra

function interval_matrices(η)
    IA = IntervalArithmetic

    # problem dimension
    n = 2 * η

    # parameters (see Table 3 in [1])
    R = IA.Interval(0.99, 1.01)
    Rdriver = IA.Interval(9.9, 10.1)
    L = IA.Interval(1e-10, 1e-10)
    C = IA.Interval(3.99e-13, 4.01e-13)

    # parameter intervals p₁-p₄
    p₁ = 1.0 / L
    p₂ = 1.0 / C
    p₃ = Rdriver / L
    p₄ = R / L

    # interval matrix that consists of four blocks A = [A₁₁ A₁₂; A₂₁ A₂₂] where
    # A₁₁ is the zero matrix,
    # A₁₂ has -p₂ on the diagonal and p₂ on the upper diagonal,
    # A₂₁ has p₁ on the diagonal and -p₁ on the lower diagonal, and
    # A₂₂ has -p₄ on the diagonal except for the top left entry, which is -p₃.
    A = A_matrix_diagonals(η, p₁, p₂, p₃, p₄)

    # interval matrix p₁r, which is a zero column matrix except for the (η + 1)
    # entry, which is p₁ (paper) resp. -p₁ (CORA implementation; see TODO below)
    r = IntervalMatrix(fill(IA.Interval(0.0), (n, 1)))
    r[η + 1, 1] = 1.0  # TODO CORA implementation says -1 but paper says +1
    B = p₁ * r

    # correction due to time scaling
    A *= 1e-9
    B *= 1e-9

    return A, B
end

function A_matrix_diagonals(η, p₁, p₂, p₃, p₄)
    A₁₁ = zeros(η, η)
    A₁₂ = Bidiagonal(fill(-p₂, η), fill(p₂, η-1), :U)
    A₂₁ = Bidiagonal(fill(p₁, η), fill(-p₁, η-1), :L)
    A₂₂ = Diagonal(vcat(-p₃, fill(-p₄, η-1)))
    A  = [A₁₁ A₁₂; A₂₁ A₂₂]
end

function A_matrix_dense(η, n, p₁, p₂, p₃, p₄)
    A = zeros(n, n)
    A[η + 1, 1] = p₁
    A[η, n] = -p₂
    A[η + 1, η + 1] = -p₃
    for i in 2:η
        A[η + i, i-1] = -p₁
        A[η + i, i] = p₁
        A[η + i, η + i] = -p₄
    end
    for i in 1:(η - 1)
        A[i, η + i] = -p₂
        A[i, η + i + 1] = p₂
    end
end

function A_matrix_paper(η, n, p₁, p₂, p₃, p₄)
    # interval matrices Q₁-Q₄
    # state dimensions: [U₁, …, U_η, I₁, …, I_η]
    Z = IntervalMatrix(fill(IntervalArithmetic.Interval(0.0), (n, n)))
    Q₁ = copy(Z)
    Q₂ = copy(Z)
    Q₃ = copy(Z)
    Q₄ = copy(Z)
    Q₁[η + 1, 1] = 1.0
    Q₂[η, n] = -1.0
    Q₃[η + 1, η + 1] = -1.0
    for i in 2:η
        Q₁[η + i, i-1] = -1.0
        Q₁[η + i, i] = 1.0
        Q₄[η + i, η + i] = -1.0
    end
    for i in 1:(η - 1)
        Q₂[i, η + i] = -1.0
        Q₂[i, η + i + 1] = 1.0
    end

    # interval matrix p₁Q₁ + p₂Q₂ + p₃Q₃ + p₄Q₄ (see Eq. (15) in [1])
    A = p₁ * Q₁ + p₂ * Q₂ + p₃ * Q₃ + p₄ * Q₄
end

function transmission_line_model(η::Int=20)
    # system and input matrices
    A, B = interval_matrices(η)

    # input domain
    U = ConstantInput(LazySets.Interval(0.99, 1.01))

    # continuous LTI system
    S = ConstrainedLinearControlContinuousSystem(A, B, nothing, U)

    return S
end
