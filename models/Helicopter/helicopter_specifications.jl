using LazySets, MathematicalPredicates, SparseArrays

function helicopter_specification()
    # initial set:
    # xᵢ ∈ [-0.1, 0.1] for i ∈ [1, 8]
    # xᵢ = 0 otherwise
    X0 = Hyperrectangle(spzeros(28), sparsevec(1:8, fill(0.1, 8), 28))

    # safety property: x1 ≤ 0.12
    property = is_contained_in(HalfSpace(sparsevec([1], [1.0], 28), 0.12))

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Dict(:T => time_horizon, :property => property)

    return X0, O
end
