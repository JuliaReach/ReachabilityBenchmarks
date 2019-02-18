include("laubloomis.jl")

using TaylorModels
#using TaylorIntegration, TaylorSeries, IntervalArithmetic

using LinearAlgebra: norm

# Equations of motion
@taylorize function laubloomis!(t, x, dx)
    dx[1] = 1.4*x[3] - 0.9*x[1]
    dx[2] = 2.5*x[5] - 1.5*x[2]
    dx[3] = 0.6*x[7] - 0.8*x[2]*x[3]
    dx[4] = 2 - 1.3*x[3]*x[4]
    dx[5] = 0.7*x[1] - x[4]*x[5]
    dx[6] = 0.3*x[1] - 3.1*x[6]
    dx[7] = 1.8*x[6] - 1.6*x[2]*x[7]
    return dx
end

# TODO: use vanderpol.jl and wrap as an algo

"""
    laubloomis_TMJets(; [t0], [T], [W], [abs_tol], [orderT], [orderQ], [maxsteps], [float_coeffs])

Build and run the Laub-Loomis model.

### Input

- `t0`       -- (optional, default: `0.0`) initial time
- `T`        -- (optional, default: `20.0`) time horizon
- `W`        -- (optional, default: `0.01`) width of the initial states
- `abs_tol`  -- (optional, default: `1e-20`) absolute tolerance used for time step
- `orderT`   -- (optional, default: `18`) order of the Taylor model in t
- `orderQ`   -- (optional, default: `9`) order of the Taylor model for Jet transport
                variables
- `maxsteps` -- (optional, default: `200`) use this maximum number of steps in
                the validated integration
- `float_coeffs` -- (optional, default: `true`) if `true`, use floating point numbers
                    for the coefficients of the polynomial variables; otherwise
                    use intervals
"""
function laubloomis_TMJets(; t0=0.0, T=20.0, W=0.01, abs_tol=1e-20, orderT=18, orderQ=9,
                            maxsteps=200, float_coeffs=true)

    # Initial conditions as mid-point of provided intervals
    q0 = IntervalBox(1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45)
    if float_coeffs
        # converts the IntervalBox into a 2-dimensional (static) array,
        # the center of the box, in this case (1.4, 2.4)
       q0 = mid.(q0)
    end

    # initial box (around `q0`) of the initial conditions
    δq0 = IntervalBox(-W..W, -W..W, -W..W, -W..W, -W..W, -W..W, -W..W)

    # TODO: wrap as a Reachability algorithm
    tTM, xTM = validated_integ(laubloomis!, q0, δq0, t0, T, orderQ, orderT, abs_tol, maxsteps=maxsteps)
    return tTM, xTM
end

# return `true` if the specification is specification is satisfied and `false`
# otherwise 
function check_property(xTM)
    satisfied = true
    for ind in eachindex(xTM[:])
        # check if solution crosses unsafe set (x4 ≥ 4.5)
        if (xTM[ind][4] >= 4.5)
            satisfied = false
            break
        end
    end
    return satisfied
end

#=
# TODO : outsource
using Plots

# plot solution and spec
function plot_solution(xTM)
    p = plot(xTM[:], color=:lightblue, legend=:false)
    plot!(p, [-2.5, 2.5], [2.75, 2.75], color=:red)
    xlims!(p, -2.5,2.5)
    ylims!(p, -3,3)
    return p
end
=#
