using LazySets, MathematicalPredicates, SparseArrays

function mna1_specification()
    # initial set: xᵢ ∈ [0.001, 0.0015] if i ≤ 2 and xᵢ = 0 otherwise
    X0 = Hyperrectangle(; low=[fill(0.001, 2); zeros(576)],
                        high=[fill(0.0015, 2); zeros(576)])

    # safety property: x1 ≤ 0.5
    property = is_contained_in(HalfSpace(sparsevec([1], [1.0], 578), 0.5))

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Dict(:T => time_horizon, :property => property)

    return X0, O
end
