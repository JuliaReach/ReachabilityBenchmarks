include("vanderpol.jl")

using TaylorIntegration, TaylorSeries, IntervalArithmetic

function vanderpol_TMJets(; k=6)
    variables = [TaylorN(1, order=2), TaylorN(2, order=2)]

    # define the set of initial states X₀ = {x: V₀(x) <= 0}
    V₀ = (x₁ - 1.4)^2 + (x₂ - 2.4)^2  - 0.15
    X0 = @set V₀ <= 0

    (𝑃, 𝑂) = vanderpol(X0=X0, variables = (x₁, x₂))

    # algorithm-specific options
    𝑂_XFZ18 = Options()
    
    # constraints Y = {x: g(x) >= 0} compact search space Y x [0, T]
    g = 25 - x₁^2 - x₂^2
    𝑂_XFZ18[:search_space] = @set g >= 0

    # degree of the relaxation
    𝑂_XFZ18[:relaxation_degree] = k

    # define the optimization solver
    𝑂_TMJets[:solver] = MosekOptimizer

    return (𝑃, 𝑂, 𝑂_TMJets)
end
