#=
Model: building.jl

See also:
bak2017cav_repeatability or
[Building example in Hylaa](https://github.com/stanleybak/hylaa/blob/master/examples/building/building.py)
=#
using Reachability, MAT, Plots

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
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

    # ===============
    # Problem solving
    # ===============

    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => LinearConstraintProperty([1., 0.], 6e-3), # x25 < 6e-3
        :vars => [25], # variable needed for property
        :partition => [(2*i-1:2*i) for i in 1:24], # 2D blocks
        :plot_vars => [0, 25]
        ), Options(input_options...))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@filename_to_png))
        toc()
    end
end # function

compute(:δ => 0.002, :N => 3); # warm-up
compute(:δ => 0.002, :T => 20.0); # benchmark settings (long)
