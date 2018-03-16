#=
Model: fom.jl

This is a 1006 x 1006 dimensional model with 1 input.
=#
using Reachability, MAT, Plots

compute(o::Pair{Symbol,<:Any}...) = compute(Options(Dict{Symbol,Any}(o)))

function compute(input_options::Options)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "fom.mat")
    A = read(file, "A")*1.0  # this model provides a matrix with Int components

    # initial set
    X0 = Hyperrectangle(zeros(1006), [fill(0.0001, 400); zeros(606)])

    # input set
    B = read(file, "B")
    U = B * BallInf([0.0], 1.0)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============
    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [1],
                                  :partition => [(2*i-1:2*i) for i in 1:503], # 2D blocks
                                  :plot_vars => [0, 1],
                                  :assume_sparse => true)
        # :projection_matrix => sparse(read(matopen(@relpath "out.mat"), "M")),
    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => 1:1006, # variables needed for property
                                  :partition => [(2*i-1:2*i) for i in 1:503], # 2D blocks
                                  :property => LinearConstraintProperty(read(matopen(@relpath "out.mat"), "M")[1,:], 185.), # y < 185
                                  :assume_sparse => true)
    end

    result = solve(S, merge(problem_options, input_options))

    # ========
    # Plotting
    # ========
    if input_options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@relpath "fom.png"))
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
        :set_type=>Interval, :partition => [[i] for i in 1:1006]);
