# ===========================================================================
# Spacecraft Rendezvous model
# See https://easychair.org/publications/paper/S2V for a description and
# https://gitlab.com/goranf/ARCH-COMP/blob/master/2018/NLN/C2E2/ARCH18_NLN/rendezvous/c2e2_nonlinear_pass_4d.hyxml
# for a reference model
# ===========================================================================
using SparseArrays, HybridSystems, MathematicalSystems, LazySets, Reachability, TaylorIntegration
using Reachability: solve

const μ = 3.986e14 * 60^2
const r = 42164e3
const mc = 500.0
const n = sqrt(μ / r^3)

const K₁ = [-28.8287 0.1005 -1449.9754 0.0046;
            -0.087 -33.2562 0.00462 -1451.5013]
const K₂ = [-288.0288 0.1312 -9614.9898 0.0;
            -0.1312 -288.0 0.0 -9614.9883]

# dynamics
@taylorize function spacecraft!(t, x, dx)
    local rc = sqrt((r + x[1])^2 + x[2]^2)
    local uxy = K₁ * x

    # x' = vx
    dx[1] = x[3]

    # y' = vy
    dx[2] = x[4]

    # vx' = n²x + 2n*vy + μ/(r^2) + μ/(rc^3)*(r+x) + ux/mc
    dx[3] = (n^2*x[1] + 2*(n*x[4])) + ((μ/(r^2) + μ/(rc^3)*(r + x[1])) + uxy[1]/mc)

    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    dx[4] = (n^2*x[2] - 2*(n*x[3])) - (μ/(rc^3)*x[2] - uxy[2]/mc)

    return dx
end

# dynamics in the aborting mode
@taylorize function spacecraft_aborting!(t, x, dx)
    local rc = sqrt((r + x[1])^2 + x[2]^2)

    # x' = vx
    dx[1] = x[3]

    # y' = vy
    dx[2] = x[4]

    # vx' = n²x + 2n*vy + μ/(r^2) + μ/(rc^3)*(r+x)
    dx[3] = (n^2*x[1] + 2*n*x[4]) + μ/(r^2) + μ/(rc^3)*(r + x[1])

    # vy' = n²y - 2n*vx - μ/(rc^3)y
    dx[4] = (n^2*x[2] - 2*n*x[3]) - μ/(rc^3)*x[2]

    return dx
end

function spacecraft_TMJets(property; T=200.0, orderT=10, orderQ=2, abs_tol=1e-10, max_steps=500)
    X0 = Hyperrectangle([-900., -400., 0., 0.], [25, 25, 0, 0.])
    𝐹 = BlackBoxContinuousSystem(spacecraft!, 4)
    𝑃 = InitialValueProblem(𝐹, X0)
    𝑂 = Options(:property=>property, :T=>T)
    𝑂jets = Options(:orderT=>orderT, :orderQ=>orderQ, :abs_tol=>abs_tol, :max_steps=>max_steps)
    return 𝑃, 𝑂, 𝑂jets
end
