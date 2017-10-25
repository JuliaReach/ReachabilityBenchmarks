#=
Model: template_model.jl

Currently, the coeffs are from the fiveDimSys.
=#
using Reachability, LazySets

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    A = sparse([-1 -4 0 0 0 0; 4 -1 0 0 0 0; 0 0 -3 1 0 0; 0 0 -1 -3 0 0; 0 0 0 0 -2. 0; 0 0 0 0 0 0])

    # initial set
    X0 = BallInf(ones(size(A, 1)), 0.1)

    # input set
    U = Hyperrectangle(zeros(6), [0.2/2, 0.5/2, 0.2/2, 0.5/2, 0.5/2, 0])

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

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
