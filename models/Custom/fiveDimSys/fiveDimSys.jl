#=
Model:  fiveDimSys.jl

This is a five-dimensional model taken from [Gir05]. This also appears as example
4.1. page 57 in the thesis [].

See also Fig. 5.1. page 71 in the same thesis.

[Gir05] -- Reachability of Linear Systems using support functions 
[ColasLeGuernicThesis] -- Reachability of Linear Systems using support functions
=#
using Reachability, LazySets

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    D = sparse([-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2.])
    P = sparse([0.6 -0.1 0.1 0.7 -0.2; -0.5 0.7 -0.1 -0.8 0; 0.9 -0.5 0.3 -0.6 0.1;
        0.5 -0.7 0.5 0.6 0.3; 0.8 0.7 0.6 -0.3 0.2])
    A = P * D * inv(full(P))
    A = sparse(A)
    A = add_spare_dimension(A)

    # initial set
    X0 = BallInf(ones(size(A, 1)), 0.1)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.01, # time step
        :blocks => [1],
        :plot_vars => [0, 1]
        ), Options(Dict{Symbol,Any}(input_options)))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        println("Plotting...")
        tic()
        options_plot = Options(
            :plot_vars => options[:plot_vars],
            :plot_name => @filename_to_png
#           :plot_indices => range_last_x_percent(length(result), 10, 3)
            )
        plot(result, options_plot)
        toc()
    end
end # function
nothing
