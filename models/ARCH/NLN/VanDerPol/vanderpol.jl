# =================================================================
# Van der Pol model
# See https://easychair.org/publications/paper/gjfh
# =================================================================

using Reachability, HybridSystems, MathematicalSystems, LazySets

function van_der_pol()

    # polynomial variables
    @polyvar x₁ x₂

    # instantiate the polynomial system
    f = [x₂, -0.2*x₁ + x₂ - 0.2*x₁^2*x₂]
    𝑃 = PolynomialContinuousSystem(f)

    # define the set of initial states X₀ = {x: V₀(x) <= 0}
    V₀ = x₁^2 + x₂^2 - 0.25
    X0 = @set V₀ <= 0

    # instantiate the IVP
    𝒮 = InitialValueProblem(𝑃, X0);

    𝑂 = Options()

    # time horizon
    𝑂[:T] = 2.0

    # variables to comptute and to plot
    𝑂[:vars] = [1, 2]
    𝑂[:plot_vars] = [1, 2]

    return (𝒮, 𝑂)
end

#=
these are algorithm-specific options

    # constraints Y = {x: g(x) >= 0} compact search space Y x [0, T]
    g = 25 - x₁^2 - x₂^2
    𝑂[:search_space] = g

    # degree of the relaxation
    k = 6
    𝑂[:relaxation_degree] = k

    # define the optimization solver
    𝑂[:solver] = MosekOptimizer

=#
