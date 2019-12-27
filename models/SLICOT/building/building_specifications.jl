using LazySets, Reachability, SparseArrays

function building_specification()
    # initial set:
    # - xᵢ ∈ [0.0002, 0.00025] for i ∈ [1, 10]
    # - xᵢ ∈ [-0.0001, 0.0001] for i = 25
    # - xᵢ = 0 for all other i
    X0 = Hyperrectangle(low=[fill(0.0002, 10); zeros(14); -0.0001; zeros(23)],
                        high=[fill(0.00025, 10); zeros(14); 0.0001; zeros(23)])

    # safety property: x25 ≤ 6e-3
    property = SafeStatesProperty(HalfSpace(sparsevec([25], [1.0], 48), 6e-3))

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Options(:T => time_horizon, :property => property)

    return X0, O
end
