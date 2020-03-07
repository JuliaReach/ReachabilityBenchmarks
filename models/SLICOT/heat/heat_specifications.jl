using LazySets, MathematicalPredicates, SparseArrays

function heat_specification()
    # initial set: xᵢ ∈ [0.6, 0.625] for i ≤ 2 and xᵢ = 0 for i > 2
    X0 = Hyperrectangle(low=[fill(0.6, 2); zeros(198)],
                        high=[fill(0.625, 2); zeros(198)])

    # safety property: x133 ≤ 0.1
    property = is_contained_in(HalfSpace(sparsevec([133], [1.0], 200), 0.1))

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Dict(:T => time_horizon, :property => property)

    return X0, O
end
