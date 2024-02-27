include("vanderpol.jl")

using DynamicPolynomials, SemialgebraicSets, MathOptInterfaceMosek

# TODO: do we return the triple (𝑃, 𝑂, 𝑂_XFZ18) or rather a problem instance 𝑃
# generic options 𝑂 and a continuous post solver instance with the given specific
# options, XFZ18(𝑂_XFZ18)

function vanderpol_XFZ18(; k=6)
    @polyvar x₁ x₂

    # define the set of initial states X₀ = {x: V₀(x) <= 0}
    V₀ = (x₁ - 1.4)^2 + (x₂ - 2.4)^2 - 0.15
    X0 = @set V₀ <= 0

    (𝑃, 𝑂) = vanderpol(; X0=X0, variables=(x₁, x₂))

    # algorithm-specific options
    𝑂_XFZ18 = Options()

    # constraints Y = {x: g(x) >= 0} compact search space Y x [0, T]
    g = 25 - x₁^2 - x₂^2
    𝑂_XFZ18[:search_space] = @set g >= 0

    # degree of the relaxation
    𝑂_XFZ18[:relaxation_degree] = k

    # define the optimization solver
    𝑂_XFZ18[:solver] = MosekOptimizer

    return (𝑃, 𝑂, 𝑂_XFZ18)
end
