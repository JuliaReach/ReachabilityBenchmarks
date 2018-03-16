#=
Model: iss.jl

This is a 270-variable model of component 1r (Russian service module) of the
International Space Station (ISS).

The corresponding SpaceEx model and configuration file are iss.xml and iss.cfg.
=#
using Reachability, MAT, Plots

compute(o::Pair{Symbol,<:Any}...) = compute(Options(Dict{Symbol,Any}(o)))

function compute(input_options::Options)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "iss.mat")
    A = sparse(read(file, "A"))

    # initial set
    X0 = BallInf(zeros(size(A, 1)), .0001) # -0.0001 <= xi <= 0.0001 for all i

    # input set
    #Uraw = BallInf([0.05], .05) * BallInf([0.9], .1) * BallInf([0.95], .05)
    #Uraw = Hyperrectangle([0.05, 0.9, 0.95], [0.05, 0.1, 0.05])
    Uraw = CartesianProductArray([BallInf([0.05], .05), BallInf([0.9], .1), BallInf([0.95], .05)])
    B = read(file, "B")
    # input U is constant
    U = B * Uraw

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # property: y < 7e-4
    p = LinearConstraintProperty(read(matopen(@relpath "out.mat"), "M")[1,:], 7e-4)

    # ===============
    # Problem solving
    # ===============
    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [182],
                                  :partition => [(2*i-1:2*i) for i in 1:135], # 2D blocks
                                  :plot_vars => [0, 182],
                                  :assume_sparse => true)
    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => 136:270, # variables needed for property
                                  :partition => [(2*i-1:2*i) for i in 1:135], # 2D blocks
                                  :property => p,
                                  :assume_sparse => true)
    end

    result = solve(S, merge(problem_options, input_options))

    # ========
    # Plotting
    # ========
    if input_options[:mode] == "reach"
        println("Plotting...")
        tic()
        #project_output = options[:projection_matrix] != nothing
        #:plot_labels => add_plot_labels(options[:plot_vars], project_output)
        plot(result)
        @eval(savefig(@relpath "iss.png"))
        toc()
    end
end # function

# ===================================
# Reach tube computation, dense time
# ===================================

info("warm-up run"; prefix=" ")
compute(:δ => 1e-3, :N => 3, :mode=>"reach", :verbosity => "warn");

info("dense time, 2D blocks Hyperrectangle"; prefix="BENCHMARK SETTINGS: ")
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach");

info("dense time, 2D blocks HPolygon, cf. Table 1 HSCC"; prefix="BENCHMARK SETTINGS: ")
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach",
        :set_type=>HPolygon, :lazy_sih=>false, :ε=>Inf);

info("dense time, 1D blocks Interval"; prefix="BENCHMARK SETTINGS: ")
compute(:δ => 1e-3, :T => 20.0, :mode=>"reach",
        :set_type=>Interval, :partition => [[i] for i in 1:270]);
