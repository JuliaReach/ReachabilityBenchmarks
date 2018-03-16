#=
Model: heat.jl

This is a 200-variable model.
=#
using Reachability, MAT, Plots

compute(o::Pair{Symbol,<:Any}...) = compute(Options(Dict{Symbol,Any}(o)))

function compute(input_options::Options)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "heat.mat")
    A = read(file, "A")

    # initial set
    # - x1-x300 are 0.0,
    # - the rest is in [0.002, 0.0015]
    X0 = Hyperrectangle([fill(0.6125, 2); zeros(198)], [fill(0.0125, 2); zeros(198)])

    # input set
    B = sparse([67], [1], [1.0], size(A, 1), 1)
    U = B * BallInf([0.0], 0.5)

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # ===============
    # Problem solving
    # ===============
    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [133],
                                  :partition => [(2*i-1:2*i) for i in 1:100], # 2D blocks
                                  :plot_vars => [0, 133])
    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => [133], # variables needed for property
                                  :partition => [(2*i-1:2*i) for i in 1:100], # 2D blocks
                                  :property => LinearConstraintProperty(sparsevec([133], [1.0], 200), 0.1)) # x133 < 0.1
    end

    result = solve(S, merge(input_options, problem_options))

    # ========
    # Plotting
    # ========
    if input_options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@relpath "heat.png"))
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
        :set_type=>Interval, :partition => [[i] for i in 1:200]);
