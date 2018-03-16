#=
Model: MNA_5.jl

This is a 10913-variable model.
=#
using Reachability, MAT, Plots

compute(o::Pair{Symbol,<:Any}...) = compute(Options(Dict{Symbol,Any}(o)))

function compute(input_options::Options)
    # =====================
    # Problem specification
    # =====================
    file = matopen(@relpath "mna5.mat")
    A = sparse(read(file, "A"))

    # initial set:
    X0 = Hyperrectangle([fill(0.000225, 10); zeros(10903)], [fill(0.000025, 10); zeros(10903)])

    # input set
    B = sparse(19:27, 1:9, fill(-1., 9), size(A, 1), 9)
    U = B * Hyperrectangle([fill(0.1, 5); fill(0.2, 4)], zeros(9))

    # instantiate continuous LTI system
    S = ContinuousSystem(A, X0, U)

    # safety property: x1 < 0.2 && x2 < 0.15
    p = LinearConstraintProperty([
        Clause([LinearConstraint(sparsevec([1], [1.0], 10913), 0.2)]),
        Clause([LinearConstraint(sparsevec([2], [1.0], 10913), 0.15)])
        ])
    # ===============
    # Problem solving
    # ===============
    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [1],
                                 :partition => vcat([(2*i-1:2*i) for i in 1:5456], [10913:10913]), # 2D blocks except last (1D)
                                  :plot_vars => [0, 1],
                                  :assume_sparse => true,
                                  :lazy_expm => true)

    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => 1:2, # variables needed for property
                                  :partition => vcat([(2*i-1:2*i) for i in 1:5456], [10913:10913]), # 2D blocks except last (1D)
                                  :property => p
                                  :assume_sparse => true,
                                  :lazy_expm => true)
    end

    result = solve(S, merge(problem_options, input_options))

    # ========
    # Plotting
    # ========
    if input_options[:mode] == "reach"
        println("Plotting...")
        tic()
        plot(result)
        @eval(savefig(@relpath "mna5.png"))
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
        :set_type=>Interval, :partition => [[i] for i in 1:10913]);
