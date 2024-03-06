using ReachabilityBenchmarks, LazySets, MathematicalPredicates, MAT

function fom_specification()
    # initial set: xᵢ ∈ ±0.0001 if i ≤ 400 and xᵢ = 0 otherwise
    X0 = Hyperrectangle(zeros(1006), [fill(0.0001, 400); zeros(606)])

    # safety property: y ≤ 185 for linear combination y (defined in out.mat)
    y = read(matopen(@current_path "out.mat"), "M")[1, :]
    property = is_contained_in(HalfSpace(y, 185.0))

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Dict(:T => time_horizon, :property => property)

    return X0, O
end
