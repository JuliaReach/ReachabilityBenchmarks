# To recreate Figure 3, run the following code from the REPL.
# Just type:
#
#     include("create_figure_3.jl")

using LazySets, Plots
if VERSION >= v"0.7"
    using LinearAlgebra
    expm = exp
else
    using Compat
end

#disable graphics output in GR - https://github.com/JuliaPlots/Plots.jl/issues/1182
ENV["GKSwstype"] = "100"

function reach_continuous(A, X0, δ, μ, T, max_order)
    # bloating factors
    Anorm = norm(A, Inf)
    α = (expm(δ*Anorm) - 1 - δ*Anorm)/norm(X0, Inf)
    β = (expm(δ*Anorm) - 1)*μ/Anorm

    # discretized system
    n = size(A, 1)
    ϕ = expm(δ*A)
    N = floor(Int, T/δ)

    # preallocate array
    R = Vector{LazySet}(undef, N)
    if N == 0
        return R
    end

    # initial reach set in the time interval [0, δ]
    ϕp = (I+ϕ)/2
    ϕm = (I-ϕ)/2
    c = X0.center
    gens = hcat(ϕp * X0.generators, ϕm * c, ϕm * X0.generators)
    R[1] = minkowski_sum(Zonotope(ϕp * c, gens),
                         Zonotope(zeros(n), Matrix((α + β)*I, n, n)))
    if order(R[1]) > max_order
        R[1] = reduce_order(R[1], max_order)
    end

    # set recurrence for [δ, 2δ], ..., [(N-1)δ, Nδ]
    ballβ = Zonotope(zeros(n), Matrix(β*I, n, n))
    for i in 2:N
        R[i] = minkowski_sum(linear_map(ϕ, R[i-1]), ballβ)
        if order(R[i]) > max_order
            R[i] = reduce_order(R[i], max_order)
        end
    end
    return R
end

function reach_hybrid(As, Ts, init, δ, μ, T, max_order)
    # start with initial states at time t=0
    queue = [(init[1], init[2], 0.)]

    res = Tuple{LazySet, Int}[]
    while !isempty(queue)
        init, loc, t = pop!(queue)
        println("currently in location $loc at time $t")
        # compute continuous successors
        R = reach_continuous(As[loc], init, δ, μ, T-t, max_order)
        found_transition = false
        for i in 1:length(R)-1
            S = R[i]
            push!(res, (S, loc))
            # check intersection with guards
            for (guard, tgt_loc) in Ts[loc]
                if !isdisjoint(S, guard)
                    # nonempty intersection with a aguard
                    new_t = t + δ * i
                    push!(queue, (S, tgt_loc, new_t))
                    found_transition = true
                    println("transition $loc -> $tgt_loc at time $new_t")
                end
            end
            # stop for first intersection
            if found_transition
                break
            end
        end
        if !found_transition && length(R) > 0
            push!(res, (R[end], loc))
        end
    end
    return res
end

function plot_res(res)
    p = plot()
    for i in 1:length(res)
        if res[i][2] == 1
            c = "blue"
        elseif res[i][2] == 2
            c = "red"
        end
        plot!(p, reduce_order(res[i][1], 2), color=c, alpha=0.1)
    end
    return p
end

function example()
    # dynamics
    A1 = [-1 -4; 4 -1]
    A2 = [1 4; -4 -1]
    As = [A1, A2]

    # transitions
    t1 = [(Hyperplane([1., 0.], -0.5), 2)]
    t2 = [(Hyperplane([0., 1.], -0.3), 1)]
    Ts = [t1, t2]

    # initial condition
    X0 = Zonotope([1.0, 0.0], Matrix(0.1*I, 2, 2))
    init_loc = 1
    init = (X0, init_loc)

    # input uncertainty
    μ = 0.001

    # discretization step
    δ = 0.001

    # time bound
    T = 4.

    # maximum order of zonotopes
    max_order = 10

    # run analysis
    res = reach_hybrid(As, Ts, init, δ, μ, T, max_order)

    # plot result
    plot_res(res)
end

example()

savefig("Figure 3.png")
