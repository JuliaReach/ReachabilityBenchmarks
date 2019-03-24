# =================================================================
# Van der Pol model
# See https://easychair.org/publications/paper/gjfh
# =================================================================

using Reachability, MathematicalSystems, LazySets, TaylorIntegration
using Reachability: solve

# Equations of motion: we write the function such that the operations are either
# unary or binary
@taylorize function vanderPol!(t, x, dx)
    local μ = 1.0
    dx[1] = x[2]
    dx[2] = (μ * x[2]) * (1 - x[1]^2) - x[1]
    return dx
end

# equations, x' = f(x(t))
F = BlackBoxContinuousSystem(vanderPol!, 2)

# initial states
X0 = Hyperrectangle(low=[1.25, 2.35], high=[1.55, 2.45])

# instantiate problem
𝑃 = InitialValueProblem(F, X0)
