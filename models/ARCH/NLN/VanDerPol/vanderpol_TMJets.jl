include("vanderpol.jl")

using TaylorModels
using TaylorModels: validated_integ

# TODO: use vanderpol.jl and wrap as an algo

"""
    vanderpol_TMJets(; [t0], [T], [abs_tol], [orderT], [orderQ],
                       [maxsteps], [property])

### Input

- `t0`       -- (optional, default: `0.0`) initial time
- `T`        -- (optional, default: `7.0`) time horizon
- `abs_tol`  -- (optional, default: `1e-10`) absolute tolerance used for time step
- `orderT`   -- (optional, default: `10`) order of the Taylor model in t
- `orderQ`   -- (optional, default: `2`) order of the Taylor model for Jet transport
                variables
- `maxsteps` -- (optional, default: `500`) use this maximum number of steps in
                the validated integration
- `property` -- (optional, default: `(t,x)->x[2] <= 2.75`) safe states property
"""
function vanderpol_TMJets(; t0=0.0, T=7.0, abs_tol=1e-10, orderT=10, orderQ=2,
                          maxsteps=500, property=(t, x) -> x[2] <= 2.75)

    # Initial conditions as mid-point of provided intervals
    q0 = [1.4, 2.4]

    # initial box (around `q0`) of the initial conditions
    δq0 = IntervalBox(-0.15 .. 0.15, -0.05 .. 0.05)

    # returns a TaylorN vector, each entry corresponding to an indep variable
    set_variables("δ"; numvars=length(q0), order=2 * orderQ)

    tTM, xTM = validated_integ(vanderPol!, q0, δq0, t0, T, orderQ, orderT,
                               abs_tol; maxsteps=maxsteps,
                               check_property=property)
    return tTM, xTM
end

#=
# return `true` if the specification is satisfied and `false`
# otherwise

function check_property(xTM)
    satisfied = true
    for ind in eachindex(xTM[:])
        # check if solution crosses unsafe set (y ≥ 2.75)
        if (xTM[ind][2] >= 2.75)
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

=#
