include("laubloomis.jl")

using TaylorModels
using TaylorModels: validated_integ

# Equations of motion
# We write the function such that the operations are either unary or binary:
@taylorize function laubloomis!(dx, x, params, t)
    dx[1] = 1.4*x[3] - 0.9*x[1]
    dx[2] = 2.5*x[5] - 1.5*x[2]
    dx[3] = 0.6*x[7] - 0.8*(x[2]*x[3])
    dx[4] = 2 - 1.3*(x[3]*x[4])
    dx[5] = 0.7*x[1] - (x[4]*x[5])
    dx[6] = 0.3*x[1] - 3.1*x[6]
    dx[7] = 1.8*x[6] - 1.6*(x[2]*x[7])
    return dx
end

"""
    laubloomis_TMJets(; [t0], [T], [W], [abs_tol], [orderT], [orderQ],
                        [maxsteps], [property])

Build and run the Laub-Loomis model.

### Input

- `t0`       -- (optional, default: `0.0`) initial time
- `T`        -- (optional, default: `20.0`) time horizon
- `W`        -- (optional, default: `0.01`) width of the initial states
- `abs_tol`  -- (optional, default: `1e-10`) absolute tolerance used for time step
- `orderT`   -- (optional, default: `7`) order of the Taylor model in t
- `orderQ`   -- (optional, default: `2`) order of the Taylor model for Jet transport
                variables
- `maxsteps` -- (optional, default: `200`) use this maximum number of steps in
                the validated integration
- `property` -- (optional, default: `(t,x) -> x[4] < 4.5`) safe states property
"""
function laubloomis_TMJets(; t0=0.0, T=20.0, W=0.01, abs_tol=1e-10,
                             orderT=7, orderQ=2, maxsteps=1000,
                             property=(t,x) -> x[4] < 4.5)

    # center of initial conditions
    q0 = [1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45]

    # initial box (around `q0`) of the initial conditions
    δq0 = IntervalBox(-W..W, 7)

    # set variables
    set_variables("δ", numvars=length(q0), order=2orderQ)

    # TODO: wrap as a Reachability algorithm
    tTM, xTM = validated_integ(laubloomis!, q0, δq0, t0, T, orderQ, orderT,
                               abs_tol, maxsteps=maxsteps,
                               check_property=property)
    return tTM, xTM
end

#=
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

using Plots

# plot solution and spec
function plot_solution(xTM)
    p = plot(xTM[:], color=:lightblue, legend=:false)
    plot!(p, [-2.5, 2.5], [2.75, 2.75], color=:red)
    xlims!(p, -2.5,2.5)
    ylims!(p, -3,3)
    return p
end

julia> data = [IntervalBox(IntervalArithmetic.Interval(tTM[i], tTM[i+1]), xTM[i][4]) for i in 1:length(tTM)-1];
julia> plot(data[:], legend=false, xlab="time", ylab="x4")

=#
