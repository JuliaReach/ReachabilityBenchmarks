# =================================================================
# Van der Pol model
# See https://easychair.org/publications/paper/gjfh
# =================================================================

using Reachability: Options, LinearConstraintProperty
using MathematicalSystems, LazySets
using DynamicPolynomials, SemialgebraicSets

"""
    vanderpol(; [T], [X0])

Construct the Van der Pol model.

### Input

- `T`  --  (optional, default: `7.0`) the time horizon for the initial
           value problem
- `X0` -- (optional, default: `[1.25, 1.55] × [2.35, 2.45]`) set of initial states

### Output

The tuple `(𝑃, 𝑂)` where `𝑃` is an initial-value problem and `𝑂` are the options.
"""
function vanderpol(; T=7.0,
                     X0=Hyperrectangle(low=[1.25, 2.35], high=[1.55, 2.45]),
                     variables=@polyvar x₁ x₂)

    𝑂 = Options()
    x₁, x₂ = variables
    𝑂[:variables] = variables
    𝑂[:vars] = [1, 2]

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
    𝑂[:property] = LinearConstraintProperty([0., 1.], 2.75)   # uses supp func evaluation
    #𝑂[:property] = SubsetProperty(HalfSpace([0., 1.], 2.75)) # uses inclusion test

    return (𝑃, 𝑂)
end
