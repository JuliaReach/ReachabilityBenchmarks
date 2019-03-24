using TaylorModels
using TaylorModels: validated_integ

# parameters of the model
const g = 9.81           # gravity constant in m/s^2
const R = 0.1            # radius of center mass in m
const l = 0.5            # distance of motors to center mass in m
const Mrotor = 0.1       # motor mass in kg
const M = 1.0            # center mass in kg
const m = M + 4*Mrotor   # total mass in kg
const mg = m*g

# moments of inertia
const Jx = (2/5)*M*R^2 + 2*l^2*Mrotor
const Jy = Jx
const Jz = (2/5)*M*R^2 + 4*l^2*Mrotor
const Cyzx = (Jy - Jz)/Jx
const Czxy = (Jz - Jx)/Jy
const Cxyz = 0.0 #(Jx - Jy)/Jz

# considering the control parameters as *parameters*
const u₁ = 1.0
const u₂ = 0.0
const u₃ = 0.0

@inline function quad_property(t, x)
    b1 = (x[3] < 1.4)
    b2 = t ≥ 1.0 ? (x[3] > 0.9) : true
    b3 = t ≥ 5.0 ? (0.98 ≤ x[3] ≤ 1.02) : true
    return b1 && b2 && b3
end

@taylorize function quadrotor!(t, x, dx)
    # unwrap the variables and the controllers; the last three are the controllers
    # x₁, x₂, x₃, x₄, x₅, x₆, x₇, x₈, x₉, x₁₀, x₁₁, x₁₂, u₁, u₂, u₃ = x
    x₁  = x[1]
    x₂  = x[2]
    x₃  = x[3]
    x₄  = x[4]
    x₅  = x[5]
    x₆  = x[6]
    x₇  = x[7]
    x₈  = x[8]
    x₉  = x[9]
    x₁₀ = x[10]
    x₁₁ = x[11]
    x₁₂ = x[12]

    # equations of the controllers
    F = (mg - 10*(x₃ - u₁)) + 3*x₆  # height control
    τϕ = -(x₇ - u₂) - x₁₀            # roll control
    τθ = -(x₈ - u₃) - x₁₁            # pitch control
    local τψ = 0.0                   # heading is uncontrolled
    #
    Tx = τϕ/Jx
    Ty = τθ/Jy
    Tz = τψ/Jz
    F_m = F/m

    # Some abbreviations
    sx7 = sin(x₇)
    cx7 = cos(x₇)
    sx8 = sin(x₈)
    cx8 = cos(x₈)
    sx9 = sin(x₉)
    cx9 = cos(x₉)
    #
    sx7sx9 = sx7*sx9
    sx7cx9 = sx7*cx9
    cx7sx9 = cx7*sx9
    cx7cx9 = cx7*cx9
    sx7cx8 = sx7*cx8
    cx7cx8 = cx7*cx8
    sx7_cx8 = sx7/cx8
    cx7_cx8 = cx7/cx8    
    #
    x4cx8 = cx8*x₄
    #
    p11 = sx7_cx8*x₁₁
    p12 = cx7_cx8*x₁₂
    xdot9 = p11 + p12

    # differential equations for the quadrotor
    #    
    dx[1] = (cx9*x4cx8 + (sx7cx9*sx8 - cx7sx9)*x₅) + (cx7cx9*sx8 + sx7sx9)*x₆
    dx[2] = (sx9*x4cx8 + (sx7sx9*sx8 + cx7cx9)*x₅) + (cx7sx9*sx8 - sx7cx9)*x₆
    dx[3] = (sx8*x₄ - sx7cx8*x₅) - cx7cx8*x₆
    dx[4] = (x₁₂*x₅ - x₁₁*x₆) - g*sx8
    dx[5] = (x₁₀*x₆ - x₁₂*x₄) + g*sx7cx8
    dx[6] = (x₁₁*x₄ - x₁₀*x₅) + (g*cx7cx8 - F_m)
    dx[7] = x₁₀ + sx8*xdot9
    dx[8] = cx7*x₁₁ - sx7*x₁₂
    dx[9] = xdot9
    dx[10] = Cyzx * (x₁₁ * x₁₂) + Tx
    dx[11] = Czxy * (x₁₀ * x₁₂) + Ty
    dx[12] = Cxyz * (x₁₀ * x₁₁) + Tz
     #
    return dx
end

"""
    quad_TMJets(; t0=0.0, T=5.0, abs_tol=1e-7, orderT=5, orderQ=2, maxsteps=500,
                  property=quad_property)

Build and run the quadrotor model.

### Input

- `t0`       -- (optional, default: `0.0`) initial time
- `T`        -- (optional, default: `5.0`) time horizon
- `abs_tol`  -- (optional, default: `1e-7`) absolute tolerance used for time step
- `orderT`   -- (optional, default: `5`) order of the Taylor model in t
- `orderQ`   -- (optional, default: `2`) order of the Taylor model for Jet transport
                variables
- `maxsteps` -- (optional, default: `1000`) use this maximum number of steps in
                the validated integration
- `property` -- (optional, default: see `quad_property`) safe states property;
                see below

### Notes

The task is to change the height from 0 [m] to 1 [m] within 5 [s]. A goal region [0.98, 1.02] of the
height x3 has to be reached within 5 [s] and the height has to stay below 1.4 for all times. After
1 [s] the height should stay above 0.9 [m]. The initial position of the quadrotor is uncertain in
all directions within [−0.4, 0.4] [m] and also the velocity is uncertain within [−0.4, 0.4] [m/s] for
all directions. All other values are initialized as 0.
"""
function quad_TMJets(; t0=0.0, T=5.0, abs_tol=1e-7,
                       orderT=5, orderQ=2, maxsteps=1000,
                       property=quad_property)

    # initial conditions, deviations and initial box
    Wpos = 0.4
    Wvel = 0.4

    q0 = zeros(12)
    δq0 = IntervalBox(-Wpos..Wpos, -Wpos..Wpos, -Wpos..Wpos, -Wvel..Wvel, -Wvel..Wvel, -Wvel..Wvel, 
                      0..0, 0..0, 0..0, 0..0, 0..0, 0..0 )

    # set variables
    set_variables("δ", numvars=length(q0), order=orderQ)

    tTM, xTM = validated_integ(quadrotor!, q0, δq0, t0, T, orderQ, orderT,
                               abs_tol, maxsteps=maxsteps,
                               check_property=property)
    return tTM, xTM
end


"""
    add_time_and_project(tTM, xTM, i)

### Input

- `tTM`  -- vector with the time slices
- `xTM`  -- vector with the interval boxes
- `i`    -- index of the dimension to be plotted, e.g. if `i = 3` then return the
            reachsets for x3 as a function of time

### Notes

For plotting, to be used as

```julia
julia> include("models/ARCH/NLN/Quadrotor/quad_TMJets.jl")

julia> tTM, xTM = quad_TMJets();

julia> x3 = add_time_and_project(tTM, xTM, 3);

julia> using Plots; pyplot()

julia> begin
            plot(x3[:], legend=:none, color=:blue, lw=0, alpha=.5, xlab="t", ylab="x3")
            ylims!(-0.5, 1.5)
            plot!([0, 5], [1.4, 1.4], color=:red, legend=:none)
            plot!([0, 5], [0.98, 0.98], color=:red, legend=:none)
            plot!([0, 5], [1.02, 1.02], color=:red, legend=:none)
       end
```
"""
function add_time_and_project(tTM, xTM, i)
    N = length(tTM) # number of reachsets
    t_x = Vector{IntervalBox{2,Float64}}(undef, N-1)
    for ind in 1:N-1
        Δt = Interval(tTM[ind], tTM[ind+1])
        Δxi = union(xTM[ind][i], xTM[ind+1][i])
        t_x[ind] = IntervalBox(Δt, Δxi)
    end
    return t_x
end
