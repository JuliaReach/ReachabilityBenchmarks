using LazySets, MathematicalPredicates, Reachability, SparseArrays

function mna5_specification()
    # initial set: xᵢ ∈ [0.0002, 0.00025] if i ≤ 10 and xᵢ = 0 otherwise
    X0 = Hyperrectangle(low=[fill(0.0002, 10); zeros(10903)],
                        high=[fill(0.00025, 10); zeros(10903)])

    # property: x1 ≤ 0.2 && x2 ≤ 0.15
    property = is_contained_in(HPolyhedron([
        HalfSpace(sparsevec([1], [1.0], 10913), 0.2),
        HalfSpace(sparsevec([2], [1.0], 10913), 0.15)]))

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Options(:T => time_horizon, :property => property)

    return X0, O
end
