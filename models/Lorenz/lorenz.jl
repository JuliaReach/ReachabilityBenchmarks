#=
Lorenz model from
https://en.wikipedia.org/wiki/Lorenz_system
=#
using Reachability, MathematicalSystems, TaylorIntegration

@taylorize function lorenz!(dx, x, params, t)
    local σ = 10.0
    local β = 8.0 / 3.0
    local ρ = 28.0
    dx[1] = σ * (x[2] - x[1])
    dx[2] = x[1] * (ρ - x[3]) - x[2]
    dx[3] = x[1] * x[2] - β * x[3]
    return dx
end

𝑆 = BlackBoxContinuousSystem(lorenz!, 3)
X0 = Hyperrectangle(low=[0.9, 0.0, 0.0], high=[1.1, 0.0, 0.0])
𝑃 = InitialValueProblem(𝑆, X0)

# reach mode
𝑂 = Options(:T=>10.0, :mode=>"reach")
sol = solve(𝑃, 𝑂, op=TMJets(:abs_tol=>1e-27, :orderT=>10, :orderQ=>2, :max_steps=>50_000));
