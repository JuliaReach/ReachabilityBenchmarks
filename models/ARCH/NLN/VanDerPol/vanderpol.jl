# =================================================================
# Van der Pol model
# See https://easychair.org/publications/paper/gjfh
# =================================================================

using Reachability, MathematicalSystems, LazySets, TaylorIntegration
using Reachability: solve

# Equations of motion: we write the function such that the operations are either
# unary or binary
@taylorize function vanderPol_mu_one!(t, x, dx)
    local μ = 1.0
    dx[1] = x[2]
    dx[2] = (μ * x[2]) * (1 - x[1]^2) - x[1]
    return dx
end

@taylorize function vanderPol_mu_two!(t, x, dx)
    local μ = 2.0
    dx[1] = x[2]
    dx[2] = (μ * x[2]) * (1 - x[1]^2) - x[1]
    return dx
end

function vanderpol(; μ=1.0,
                     X0=Hyperrectangle(low=[1.25, 2.35], high=[1.55, 2.45]),
                     T=7.0,
                     property=(t, x) -> x[2] < 2.75)

    if μ == 1.0
        f = vanderPol_mu_one!
    elseif μ == 2.0
        f = vanderPol_mu_two!
    else
        error("the value of μ = $μ is not implemented")
    end
    # equations, x' = f(x(t))
    F = BlackBoxContinuousSystem(f, 2)

    # problem options
    𝑂 = Options(:T=>T, :mode=>"check", :property=>property)
    
    # instantiate problem
    𝑃 = InitialValueProblem(F, X0)

    return 𝑃, 𝑂
end

function splitX0(H::AbstractHyperrectangle, n::Int=2, m::Int=2)
     @assert dim(H) == 2
     r = copy(radius_hyperrectangle(H))
     r[1] /= n
     r[2] /= m
     d = 2*r
     result = Vector{Hyperrectangle}()
     sizehint!(result, n*m)
     c0 = low(H) - r
     c = copy(c0)
     for i in 1:n
         c[1] += d[1]
         c[2] = c0[2]
         for j in 1:m
             c[2] += d[2]
             push!(result, Hyperrectangle(c, r))
         end
     end
     return result
end
