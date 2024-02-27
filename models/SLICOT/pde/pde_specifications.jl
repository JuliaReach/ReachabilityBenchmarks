using ReachabilityBenchmarks, LazySets, MathematicalPredicates, MAT

function pde_specification()
    # initial set:
    # xᵢ = 0 for i ≤ 64
    # xᵢ ∈ [0.001, 0.0015] for i ∈ [65, 80]
    # xᵢ ∈ [-0.002, -0.0015] for i > 80
    X0 = Hyperrectangle(; low=[zeros(64); fill(0.001, 16); fill(-0.002, 4)],
                        high=[zeros(64); fill(0.0015, 16); fill(-0.0015, 4)])

    # safety property: y ≤ 12 for linear combination y (defined in out.mat)
    y = read(matopen(@relpath "out.mat"), "M")[1, :]
    property = is_contained_in(HalfSpace(y, 12.0))

    # time horizon: 20 time units
    time_horizon = 20.0

    # specification
    O = Dict(:T => time_horizon, :property => property)

    return X0, O
end
