# =================================================================
# Laub-Loomis model
# See https://easychair.org/publications/paper/gjfh
# =================================================================

using Reachability: Options, LinearConstraintProperty
using MathematicalSystems, LazySets
using DynamicPolynomials, SemialgebraicSets

"""
    laubloomis(; [T], [X0], [W])

Construct the Laub-Loomis model.

### Input

- `T`  --  (optional, default: `20.0`) the time horizon for the initial-value
            problem
- `W`  --  (optional, default: `0.01`) width of the initial states
- `X0` --  (optional, default: an axis-aligned box centered at
           `(1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45)` and radius `W`) set of initial states
- `unsafe_bound` -- (optional, default: `4.5`) bound for the variable `x₄` that
                    defines the unsafe region, of the form `x₄ ≥ unsafe_bound`
- `variables`    -- (optional, default: `PolyVar`) the set of variables used to
                    describe the polynomial ODE

### Output

The tuple `(𝑃, 𝑂)` where `𝑃` is an initial-value problem and `𝑂` are the options.
"""
function laubloomis(; T=20.0,
                      W=0.01,
                      X0=BallInf([1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45], W,),
                      unsafe_bound=4.5,
                      variables=@polyvar x₁ x₂ x₃ x₄ x₅ x₆ x₇)

    𝑂 = Options()

    # unrwap the variables
    x₁, x₂, x₃, x₄, x₅, x₆, x₇ = variables
    𝑂[:variables] = variables
    𝑂[:vars] = [1:7;]

    # instantiate the polynomial system
    f = [1.4x₃ - 0.9x₁,
         2.5x₅ - 1.5x₂,
         0.6x₇ - 0.8x₂*x₃,
         2 - 1.3x₃*x₄,
         0.7x₁ - x₄*x₅,
         0.3x₁ - 3.1x₆,
         1.8x₆ - 1.6x₂*x₇]

    𝐹 = PolynomialContinuousSystem(f)

    # instantiate the IVP
    𝑃 = InitialValueProblem(𝐹, X0)

    # time horizon
    𝑂[:T] = T

    # variables to plot
    𝑂[:plot_vars] = [0, 4]

    # safety property
    𝑂[:property] = LinearConstraintProperty([0, 0, 0, -1., 0, 0, 0], -unsafe_bound)
    # @set x₄ ≥ 0.01, vars=(x₁, x₂, x₃, x₄, x₅, x₆, x₇)
    return (𝑃, 𝑂)
end
