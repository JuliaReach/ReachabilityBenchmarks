function get_coeffs(flow, n, m, state_vars, input_vars)

    A = Matrix{Float64}(undef, n, n)
    B = Matrix{Float64}(undef, n, m)
    c = zeros(Float64, n)

    for (i, fi) in enumerate(flow)
        #println(fi)
        RHS = convert(Basic, fi.args[2])
        #println(RHS)
        #println()
        const_term = subs(RHS, [xi=>0.0 for xi in state_vars]...,
                               [ui=>0.0 for ui in input_vars]...)

        if const_term != 0.0
            C[i] = const_term
            isaffine = true
        end
        # terms linear in the *state* variables
        A[i, :] = convert.(Float64, diff.(RHS, state_vars))

        # terms linear in the *input* variables
        B[i, :] = convert.(Float64, diff.(RHS, input_vars))
    end
    return A, B, c
end
