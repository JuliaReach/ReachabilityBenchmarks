include("quadrotor.jl")

using TaylorModels
using TaylorModels: validated_integ

# =====================
# Equations of motion
# =====================
# parameters of the model
const g = 9.81           # gravity constant in m/s^2
const R = 0.1            # radius of center mass in m
const l = 0.5            # distance of motors to center mass in m
const Mrotor = 0.1       # motor mass in kg
const M = 1.0            # center mass in kg
const m = M + 4*Mrotor   # total mass in kg

# moments of inertia
const Jx = 2/5*M*R^2 + 2*l^2*Mrotor
const Jy = Jx
const Jz = 2/5*M*R^2 + 4*l^2*Mrotor

# inputs: cannot handle non-deterministic inputs)
# we take them to be the center of inputs
const u₁ = 1.0
const u₂ = 0.0
const u₃ = 0.0

# We write the function such that the operations are either unary or binary:
@taylorize function quadrotor!(t, x, dx)

    # variables
    x₁, x₂, x₃, x₄    = x[1], x[2], x[3], x[4]
    x₅, x₆, x₇, x₈    = x[5], x[6], x[7], x[8]
    x₉, x₁₀, x₁₁, x₁₂ = x[9], x[10], x[11], x[12]

    # derivatives
    ẋ₁, ẋ₂, ẋ₃, ẋ₄    = dx[1], dx[2], dx[3], dx[4]
    ẋ₅, ẋ₆, ẋ₇, ẋ₈    = dx[5], dx[6], dx[7], dx[8]
    ẋ₉, ẋ₁₀, ẋ₁₁, ẋ₁₂ = dx[9], dx[10], dx[11], dx[12]

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

    return dx
end

"""
    quadrotor_TMJets(; [t0], [T], [W], [abs_tol], [orderT], [orderQ],
                       [maxsteps], [property])

Build and run the Laub-Loomis model.

### Input

- `t0`       -- (optional, default: `0.0`) initial time
- `T`        -- (optional, default: `5.0`) time horizon
- `abs_tol`  -- (optional, default: `1e-20`) absolute tolerance used for time step
- `orderT`   -- (optional, default: `9`) order of the Taylor model in t
- `orderQ`   -- (optional, default: `2`) order of the Taylor model for Jet transport
                variables
- `maxsteps` -- (optional, default: `200`) use this maximum number of steps in
                the validated integration
- `property` -- (optional, default: `x -> true`) safe states property
"""
function quadrotor_TMJets(; t0=0.0, T=5.0, abs_tol=1e-20,
                            orderT=9, orderQ=2, maxsteps=1000,
                            property = x -> true)

    # center of initial conditions
    q0 = zeros(6)

    # initial box (around `q0`) of the initial conditions
    δq0 = IntervalBox(-0.4..0.4, -0.4..0.4, -0.4..0.4, -0.4..0.4, -0.4..0.4, -0.4..0.4,
                       0.0..0.0, 0.0..0.0, 0.0..0.0, 0.0..0.0, 0.0..0.0, 0.0..0.0)

    # set variables
    set_variables("δ", numvars=length(q0), order=2orderQ)

    # TODO: wrap as a Reachability algorithm
    tTM, xTM = validated_integ(quadrotor!, q0, δq0, t0, T, orderQ, orderT,
                               abs_tol, maxsteps=maxsteps,
                               check_property=property)
    return tTM, xTM
end
