using Reachability, MathematicalSystems, LazySets, TaylorIntegration
using Reachability: solve

# Paremeters
const μ = 3.986e14 * 60^2
const r = 42164.0e3
const r2 = r^2
const mc = 500.0
const n2 = μ / r^3
const n = sqrt(n2)

const two_n = 2 * n
const μ_r2 = μ / r2

const K₁ = [-28.8287 0.1005 -1449.9754 0.0046;
            -0.087 -33.2562 0.00462 -1451.5013]
const K₂ = [-288.0288 0.1312 -9614.9898 0.0;
            -0.1312 -288.0 0.0 -9614.9883]

const K₁mc = K₁ / mc
const K₂mc = K₂ / mc

const zK = zero(K₁)

function switch_controller!(t, x)
    bool_time = constant_term(t) < 120.0
    x_coord = constant_term(constant_term(x[1]))
    bool_approach = -100.0 ≥ x_coord ≥ -1000.0

    if bool_time

        # Approaching mode
        if bool_approach
            return K₁mc * x
        end

        # Rendezvous attempt
        return K₂mc * x
    end

    # Aborting
    return zK * x
end

# dynamics in the 'approaching' mode
@taylorize function space_rendezvous!(dx, x, params, t)
    x_1 = r + x[1]
    x_12 = x_1^2
    x_22 = x[2]^2
    rc = sqrt(x_12 + x_22)
    rc3 = rc^3
    μ_rc3 = μ / rc3

    uxy = switch_controller!(t, x)

    # x' = vx
    dx[1] = x[3]
    # y' = vy
    dx[2] = x[4]
    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    dx[3] = (n2 * x[1] + two_n * x[4]) + ((μ_r2 - μ_rc3 * x_1) + uxy[1])
    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    dx[4] = (n2 * x[2] - two_n * x[3]) - (μ_rc3 * x[2] - uxy[2])
    # t' = 1
    # dx[5] = one(t)

    return dx
end

𝑂 = Options(:T => 200.0, :plot_vars => [1, 2], :project_reachset => true, :mode => "reach")
X0 = Hyperrectangle([-900.0, -400.0, 0.0, 0.0], [25.0, 25.0, 0.0, 0.0])
𝑃 = IVP(BlackBoxContinuousSystem(space_rendezvous!, 4), X0)
𝑂jets = Options(:orderT => 10, :orderQ => 2, :abs_tol => 1e-28, :max_steps => 11000)

# solve (long time: 569.110821 seconds (7.69 G allocations: 475.990 GiB, 20.96% gc time)
# @time sol = solve(𝑃, 𝑂, op=TMJets(𝑂jets))
