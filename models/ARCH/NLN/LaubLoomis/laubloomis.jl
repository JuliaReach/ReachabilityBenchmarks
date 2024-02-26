# =================================================================
# Laub-Loomis model
# See https://easychair.org/publications/paper/gjfh
# =================================================================

using Reachability, MathematicalSystems, TaylorIntegration

# Equations of motion
# We write the function such that the operations are either unary or binary:
@taylorize function laubloomis!(dx, x, params, t)
    dx[1] = 1.4 * x[3] - 0.9 * x[1]
    dx[2] = 2.5 * x[5] - 1.5 * x[2]
    dx[3] = 0.6 * x[7] - 0.8 * (x[2] * x[3])
    dx[4] = 2 - 1.3 * (x[3] * x[4])
    dx[5] = 0.7 * x[1] - (x[4] * x[5])
    dx[6] = 0.3 * x[1] - 3.1 * x[6]
    dx[7] = 1.8 * x[6] - 1.6 * (x[2] * x[7])
    return dx
end

function laubloomis(; T=20.0, W=0.01, plot_vars=[0, 4],
                    property=(t, x) -> x[4] < 4.5,
                    project_reachset=true)

    # equations, x' = f(x(t))
    𝐹 = BlackBoxContinuousSystem(laubloomis!, 7)

    X0c = [1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45]
    X0 = Hyperrectangle(X0c, fill(W, 7))

    # instantiate the IVP
    𝑃 = InitialValueProblem(𝐹, X0)

    # general options
    𝑂 = Options(:T => T, :plot_vars => plot_vars, :property => property,
                :project_reachset => project_reachset, :mode => "check")

    return (𝑃, 𝑂)
end
