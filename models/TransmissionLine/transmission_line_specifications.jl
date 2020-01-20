using LazySets, Reachability, MathematicalSystems, IntervalMatrices
using LinearAlgebra: Diagonal

function transmission_line_specification(S::ConstrainedLinearControlContinuousSystem)
    IA = IntervalArithmetic
    # initial set: -A⁻¹bu + □(1e-3)
    # where A and b are the midpoints of S.A and S.B, respectively

    A⁻¹ = inv(mid(S.A))
    b = S.B  # TODO paper says midpoint but CORA implementation uses range
             #      we currently cannot do matrix * mid(S.B)
    n = size(b, 1)
    u = IA.Interval(-0.2, 0.2)
    □(r) = IntervalMatrix(fill(IA.Interval(-r, r), (n, 1)))
    X0_itvmat = IntervalMatrix(-A⁻¹ * b) * u + □(1e-3)
    # convert interval matrix to zonotope
    X0 = Zonotope(zeros(n), Diagonal(vec([abs(x.hi) for x in X0_itvmat])))

    # time horizon: 0.7 time units
    time_horizon = 0.7

    # specification
    O = Options(:T => time_horizon)

    return X0, O
end
