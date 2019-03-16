#=
Model: Motor (8 variables, 2 inputs)
=#
using Reachability, SparseArrays

motor(o::Pair{Symbol, <:Any}...) = motor(Options(Dict{Symbol, Any}(o)))

function motor(input_options::Options)
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
    U = Hyperrectangle([0.23, 0.3], [0.07, 0.1])

    # instantiate continuous LTI system
    S = InitialValueProblem(
        ConstrainedLinearControlContinuousSystem(A, B, nothing, U), X0)

    # property: x1 < 0.35 || x5 < 0.45
    property = Disjunction([
        SafeStatesProperty(HalfSpace([1.; zeros(7)], 0.35)),
        SafeStatesProperty(HalfSpace([zeros(4); 1.; zeros(3)], 0.45))])

    # =======================
    # Problem default options
    # =======================
    partition = [(2*i-1:2*i) for i in 1:4] # 2D blocks

    if input_options[:mode] == "reach"
        problem_options = Options(:vars => [5],
                                  :partition => partition,
                                  :plot_vars => [0, 5])
    elseif input_options[:mode] == "check"
        problem_options = Options(:vars => [1, 5],
                                  :partition => partition,
                                  :property => property)
    end

    return (S, merge(problem_options, input_options))
end
