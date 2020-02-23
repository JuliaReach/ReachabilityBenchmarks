using LazySets, MathematicalPredicates, Reachability, SparseArrays

function cdplayer_specification()
    # initial set: xᵢ ∈ [-1, 1] for all i
    X0 = BallInf(zeros(120), 1.0)

    # safety property: 2·x₁ - 3·x₂ ≤ 450.8
    property =
        is_contained_in(HalfSpace(sparsevec([1, 2], [2.0, -3.0], 120), 450.8))

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Options(:T => time_horizon, :property => property)

    return X0, O
end
