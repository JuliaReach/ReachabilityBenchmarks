# =================================================================
# Quadrotor model
# See https://easychair.org/publications/paper/gjfh
#
#
# Physical variables:
#  x₁ : interitial (north) position
#  x₂ : intertial (east) position
#  x₃ : altitude
#  x₄ : longitudinal velocity
#  x₅ : lateral velocity
#  x₆ : vertical velocity
#  x₇ : roll angle
#  x₈ : pitch angle
#  x₉ : yaw angle
#  x₁₀ : roll rate
#  x₁₁ : pitch rate
#  x₁₂ : yaw rate
# =================================================================

using Reachability: Options, LinearConstraintProperty
using MathematicalSystems, LazySets
using DynamicPolynomials, SemialgebraicSets

"""
    quadrotor(; [T], [X0], [W])

Construct the Quadrotor model.

### Input

- `T`  --  (optional, default: `20.0`) the time horizon for the initial
           value problem
- `X0` --  (optional, default: an axis-aligned box centered at
           `(1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45)` and radius `W`) set of initial states
- `W`  --  (optional, default: `0.01`) width of the initial states
- `unsafe_bound` -- (optional, default: `4.5`) bound for the variable `x₄` that
                    defines the unsafe region, of the form `x₄ ≥ unsafe_bound`

### Output

The tuple `(𝑃, 𝑂)` where `𝑃` is an initial-value problem and `𝑂` are the options.
"""
function quadrotor(; T=20.0,
                     X0=nothing,
                     W=nothing,
                     unsafe_bound=4.5,
                     variables=@polyvar x[1:12],
                     inputs=@polyvar u₁ u₂ u₃)

    # parameters of the model
    g = 9.81           # gravity constant in m/s^2
    R = 0.1            # radius of center mass in m
    l = 0.5            # distance of motors to center mass in m
    Mrotor = 0.1       # motor mass in kg
    M = 1.             # center mass in kg
    m = M + 4*Mrotor   # total mass in kg

    # moments of inertia
    Jx = 2/5*M*R^2 + 2*l^2*Mrotor
    Jy = Jx
    Jz = 2/5*M*R^2 + 4*l^2*Mrotor

    𝑂 = Options()

    # unrwap the variables and the inputs
    x₁, x₂, x₃, x₄, x₅, x₆, x₇, x₈, x₉, x₁₀, x₁₁, x₁₂ = variables[1] # or variables
    u₁, u₂, u₃ = inputs[1]

    𝑂[:variables] = variables
    𝑂[:vars] = [1:12;]

    # equations of the controllers
    F = m*g - 10*(x₃ - u₁) + 3*x₆  # height control
    τϕ = -(x₇-u₂) - x₁₀            # roll control
    τθ = -(x₈ - u₃) - x₁₁          # pitch control
    τψ = 0                         # heading is uncontrolled

    # differential equations for the quadrotor
    ẋ₁ = cos(x₈)*cos(x₉)*x₄ + (sin(x₇)*sin(x₈)*cos(x₉) - cos(x₇)*sin(x₉))*x₅
         + (cos(x₇)*sin(x₈)*cos(x₉) + sin(x₇)*sin(x₉))*x₆
    ẋ₂ = cos(x₈)*sin(x₉)*x₄ + (sin(x₇)*sin(x₈)*sin(x₉) + cos(x₇)*cos(x₉))*x₅
         + (cos(x₇)*sin(x₈)*sin(x₉) - sin(x₇)*cos(x₉))*x₆
    ẋ₃ = sin(x₈)*x₄ - sin(x₇)*cos(x₈)*x₅ - cos(x₇)*cos(x₈)*x₆
    ẋ₄ = x₁₂*x₅ - x₁₁*x₆ - g*sin(x₈)
    ẋ₅ = x₁₀*x₆ - x₁₂*x₄ + g*cos(x₈)*sin(x₇)
    ẋ₆ = x₁₁*x₄ - x₁₀*x₅ + g*cos(x₈)*cos(x₇) - F/m
    ẋ₇ = x₁₀ + sin(x₇)*tan(x₈)*x₁₁ + cos(x₇)*tan(x₈)*x₁₂
    ẋ₈ = cos(x₇)*x₁₁ - sin(x₇)*x₁₂
    ẋ₉ = (sin(x₇)/cos(x₈))*x₁₁ + (cos(x₇)/cos(x₈))*x₁₂
    ẋ₁₀ = (Jy - Jz)/Jx * x₁₁ * x₁₂ + τϕ/Jx
    ẋ₁₁ = (Jz - Jx)/Jy * x₁₀ * x₁₂ + τθ/Jy
    ẋ₁₂ = (Jx - Jy)/Jz * x₁₀ * x₁₁ + τψ/Jz

    𝐹 = PolynomialContinuousSystem(f)

    # set default of initial states
    if X0 == nothing
        if W == nothing
            W = 0.01
        end
        X0 = Hyperrectangle([1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45], fill(W, 7))
    end

    # instantiate the IVP
    𝑃 = InitialValueProblem(𝐹, X0)

    # time horizon
    𝑂[:T] = T

    # variables to plot
    𝑂[:plot_vars] = [1, 2]

    # safety property
    𝑂[:property] = LinearConstraintProperty([0, 0, 0, -1., 0, 0, 0], -unsafe_bound)
    # @set x₄ ≥ 0.01, vars=(x₁, x₂, x₃, x₄, x₅, x₆, x₇)
    return (𝑃, 𝑂)
end
