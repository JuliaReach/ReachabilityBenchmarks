using MathematicalSystems, LazySets, MathematicalPredicates,
      DynamicPolynomials, SemialgebraicSets
using Reachability: Options

"""
    vanderpol_MP(; [T], [X0], [variables])

Construct the Van der Pol model using MultivariatePolynomials.

### Input

- `T`         -- (optional, default: `7.0`) the time horizon for the initial
                  value problem
- `X0`        -- (optional, default: `[1.25, 1.55] × [2.35, 2.45]`) set of initial states
- `variables` -- (optional, default: `PolyVar`) the set of polynomial variables that
                  are used in the equations 
### Output

The tuple `(𝑃, 𝑂)` where `𝑃` is an initial-value problem and `𝑂` are the options.
"""
function vanderpol_MP(; T=7.0,
                      X0=Hyperrectangle(; low=[1.25, 2.35], high=[1.55, 2.45]),
                      variables=@polyvar x₁ x₂)
    𝑂 = Options()
    x₁, x₂ = variables
    𝑂[:variables] = variables

    # instantiate the polynomial system
    f = [x₂, x₂ - x₁ - x₁^2 * x₂]
    𝐹 = PolynomialContinuousSystem(f)

    # instantiate the IVP
    𝑃 = InitialValueProblem(𝐹, X0)

    # time horizon
    𝑂[:T] = T

    # variables to plot
    𝑂[:plot_vars] = [1, 2]

    # safety property
    𝑂[:property] = is_contained_in(HalfSpace([0.0, 1.0], 2.75))

    return (𝑃, 𝑂)
end
