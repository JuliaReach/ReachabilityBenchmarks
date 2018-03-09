#=
Model: motor.jl

For the SpaceEx model, see [1] Motor.xml and the configuration file Motor.cfg.
Output variables are x1 and x2.

[1] BakDuggirala2017cavrepeatability/arch_benchmarks/Motor/Motor.xml
=#
using Reachability, Plots

function compute(input_options::Pair{Symbol,<:Any}...)
    # =====================
    # Problem specification
    # =====================
    I = [1, 2, 2, 3, 3, 3, 3, 4, 5, 6, 6, 7, 7, 7, 7, 8]
    J = [2, 3, 2, 1, 2, 3, 4, 1, 6, 7, 6, 5, 6, 7, 8, 5]
    vals = [1, 8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0, 1.0,
           8487.2, -1.0865, -2592.1, -21.119, -698.91, -141399.0, 1.0]
    A = sparse(I, J, vals)

    # initial set
    X0 = Hyperrectangle([0.00225, 0.0, 0.0, 0.0, 0.00125, 0.0, 0.0, 0.0],
                        [0.00025, 0.0, 0.0, 0.0, 0.00025, 0.0, 0.0, 0.0])

    # input set
    B = sparse([4, 8], [1, 2], [-1.0, -1.0])
    U = B * Hyperrectangle([0.23, 0.3], [0.07, 0.1])

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============
    # define solver-specific options
    options = merge(Options(
        :mode => "reach",
        :property => LinearConstraintProperty(
                         Clause([LinearConstraint([1.; zeros(7)], 0.35),
                                 LinearConstraint([zeros(4); 1.; zeros(3)], 0.45)])),
                                 # x1 < 0.35 || x5 < 0.45
#       :vars => [5], # variable for single block analysis
        :vars => [1, 5], # variables needed for property
        :partition=> [(2*i-1:2*i) for i in 1:4], # 2D blocks
        :plot_vars => [0, 5]
        ), Options(input_options...))

    result = solve(S, options)

    # ========
    # Plotting
    # ========
    if options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig("motor.jl"))
        toc()
    end
end # function

compute(:δ => 0.001, :N => 3); # warm-up
compute(:δ => 0.001, :T => 20.0); # benchmark settings (long)
