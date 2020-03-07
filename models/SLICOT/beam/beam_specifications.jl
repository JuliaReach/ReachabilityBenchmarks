using LazySets, MathematicalPredicates, SparseArrays

function beam_specification()
    # initial set: xᵢ = 0 if i ≤ 300 and xᵢ ∈ [0.0015, 0.002] otherwise
    X0 = Hyperrectangle(low=[zeros(300); fill(0.0015, 48)],
                        high=[zeros(300); fill(0.002, 48)])

    # safety property: x89 ≤ 2100
    property = is_contained_in(HalfSpace(sparsevec([89], [1.0], 348), 2100.0))

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Dict(:T => time_horizon, :property => property)

    return X0, O
end
