using LazySets, Reachability, MathematicalSystems, IntervalMatrices

function transmission_line_specification(S::ConstrainedLinearControlContinuousSystem)
    # initial set: -A⁻¹bu + □(1e-3)
    # where A and b are the midpoints of S.A and S.B, respectively

    # TODO CORA first creates an interval matrix and then converts to zonotope.
    # is this equivalent?
    A⁻¹ = inv(mid(S.A))
    b = mid(S.B)  # TODO paper says midpoint but CORA implementation uses range
    n = size(b, 1)
    u = LazySets.Interval(-0.2, 0.2)
    □(r) = BallInf(zeros(n), r)
    X0 = minkowski_sum(linear_map(-A⁻¹ * b, u), □(1e-3))

    # time horizon: 0.7 time units
    time_horizon = 0.7

    # specification
    O = Options(:T => time_horizon)

    return X0, O
end
