include("vanderpol.jl")

using TaylorModels
using TaylorModels: validated_integ, validated_integ2

@taylorize function vanderPol!(t, x, dx)
    local μ = 1.0
    dx[1] = x[2]
    dx[2] = (μ * x[2]) * (1 - x[1]^2) - x[1]
    return dx
end

# TODO: use vanderpol.jl and wrap as an algo

"""
    vanderpol_TMJets(; [t0], [T], [abs_tol], [orderT], [orderQ], [maxsteps], [float_coeffs])

### Input

- `t0`       -- (optional, default: `0.0`) initial time
- `T`        -- (optional, default: `7.0`) time horizon
- `abs_tol`  -- (optional, default: `1e-20`) absolute tolerance used for time step
- `orderT`   -- (optional, default: `13`) order of the Taylor model in t
- `orderQ`   -- (optional, default: `9`) order of the Taylor model for Jet transport
                variables
- `maxsteps` -- (optional, default: `500`) use this maximum number of steps in
                the validated integration
- `float_coeffs` -- (optional, default: `true`) if `true`, use floating point numbers
                    for the coefficients of the polynomial variables; otherwise
                    use intervals
- `sym_norm`
- `check_property`
- `jet_transport`
"""
function vanderpol_TMJets(; t0=0.0, T=7.0, abs_tol=1e-20, orderT=13, orderQ=9,
                            maxsteps=500, float_coeffs=true, sym_norm=true,
                            check_property=x->x[2] <= 2.75,
                            jet_transport=true)

    # Initial conditions as mid-point of provided intervals
    q0 = IntervalBox(1.4, 2.4)
    if float_coeffs
         # converts the IntervalBox into a 2-dimensional (static) array,
         # the center of the box, in this case (1.4, 2.4)
        q0 = mid.(q0)
    end

    # initial box (around `q0`) of the initial conditions
    δq0 = IntervalBox(-0.15..0.15, -0.05..0.05)

    # returns a TaylorN vector, each entry corresponding to an indep variable
    set_variables("x", numvars=2, order=2*orderQ)

    # TODO: wrap as a Reachability algorithm
    if jet_transport
        tTM, xTM = validated_integ2(vanderPol!, q0, δq0, t0, T, orderQ, orderT,
                                   abs_tol, maxsteps=maxsteps,
                                   check_property=check_property,
                                   sym_norm=sym_norm)
    else
        tTM, xTM = validated_integ(vanderPol!, q0, δq0, t0, T, orderQ, orderT,
                                   abs_tol, maxsteps=maxsteps,
                                   check_property=check_property,
                                   sym_norm=sym_norm)
    end
    return tTM, xTM
end

# return `true` if the specification is specification is satisfied and `false`
# otherwise 
#=
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
=#

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
