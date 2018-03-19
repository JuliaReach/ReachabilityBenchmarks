#=
Model: motor.jl

For the SpaceEx model, see [1] Motor.xml and the configuration file Motor.cfg.
Output variables are x1 and x2.

[1] BakDuggirala2017cavrepeatability/arch_benchmarks/Motor/Motor.xml
=#
using Reachability, Plots

compute(o::Pair{Symbol,<:Any}...) = compute(Options(Dict{Symbol,Any}(o)))

function compute(input_options::Options)
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

    # property: x1 < 0.35 || x5 < 0.45
    p = LinearConstraintProperty(Clause([LinearConstraint([1.; zeros(7)], 0.35),
                                         LinearConstraint([zeros(4); 1.; zeros(3)], 0.45)]))

    # ===============
    # Problem solving
    # ===============
    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [5],
                                  :partition => [(2*i-1:2*i) for i in 1:4], # 2D blocks
                                  :plot_vars => [0, 5])
    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => [1, 5], # variables needed for property
                                  :partition => [(2*i-1:2*i) for i in 1:4], # 2D blocks
                                  :property => p)
    end

    result = solve(S, merge(problem_options, input_options))

    # ========
    # Plotting
    # ========
    if input_options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@relpath "motor.png"))
        toc()
    end
end # function

# ===================================
# Reach tube computation, dense time
# ===================================

info("warm-up run"; prefix=" ")
compute(:δ => 1e-3, :N => 3, :mode=>"reach", :verbosity => "warn");

info("dense time, 2D blocks Hyperrectangle"; prefix="BENCHMARK SETTINGS: ")
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach", :verbosity => "info");

info("dense time, 2D blocks HPolygon, cf. Table 1 HSCC"; prefix="BENCHMARK SETTINGS: ")
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach", :verbosity => "info",
        :set_type=>HPolygon, :lazy_sih=>false, :ε=>Inf);

info("dense time, 1D blocks Interval"; prefix="BENCHMARK SETTINGS: ")
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach", :verbosity => "info",
        :set_type=>Interval, :partition => [[i] for i in 1:8]);
