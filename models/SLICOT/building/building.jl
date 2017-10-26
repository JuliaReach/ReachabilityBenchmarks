#=
Model: building.jl

See also:
bak2017cav_repeatability or
[Building example in Hylaa](https://github.com/stanleybak/hylaa/blob/master/examples/building/building.py)
=#
using Reachability, LazySets, MAT

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    println("System construction...")
    tic()

    file = matopen(@relpath "building.mat")
    A = sparse(read(file, "A"))

    # initial set
    # - x1-x10 are in [0.0002, 0.00025],
    # - x25 is in [-0.0001, 0.0001],
    # - the rest is 0.0
    # 1) Hyperrectangle
    X0 = Hyperrectangle([fill(0.000225, 10); zeros(38)], [fill(0.000025, 10); zeros(14); 0.0001; zeros(23)])
    # 2) Cartesian product of intervals
    # X0 = CartesianProductArray([BallInf(fill(0.000225, 10), 0.000025), Singleton(zeros(14)), BallInf([0.], 0.0001), Singleton(zeros(23))])

    # input set
    B = sparse(read(file, "B"))
    U = B * BallInf([0.9], .1)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    toc()

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => Property([1., 0.], 6e-3), # x25 < 6e-3
        :T => 20., # time horizon
        :N => 3, # number of time steps
#       :δ => 0.002, # time step
        :blocks => [@block_id(25)],
        :plot_vars => [0, 25]
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
