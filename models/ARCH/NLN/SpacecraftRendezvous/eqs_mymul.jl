const μ = 3.986e14 * 60^2
const r = 42164.0e3
const r2 = r^2
const mc = 500.0
const n2 = μ / r^3
const n = sqrt(n2)

const two_n = 2*n
const μ_r2 = μ/r2

const K₁ = [-28.8287 0.1005 -1449.9754 0.0046;
            -0.087 -33.2562 0.00462 -1451.5013]
const K₂ = [-288.0288 0.1312 -9614.9898 0.0;
            -0.1312 -288.0 0.0 -9614.9883]

const K₁mc = K₁/mc
const K₂mc = K₂/mc

function mymul!(v, A, x)
    for ind = 1:length(v)
        v[ind] = zero(x[1])
        for jind = 1:length(x)
            v[ind] += A[ind, jind] * x[jind]
        end
    end
    return nothing
end

# dynamics in the 'approaching' mode
@taylorize function spacecraft_approaching!(t, x, dx)
    x_1 = r + x[1]
    x_12 = x_1^2
    x_22 = x[2]^2
    rc = sqrt(x_12 + x_22)
    rc3 = rc^3
    μ_rc3 = μ / rc3

    uxy = Array{typeof(x[1])}(undef, size(K₁mc,1))
    mymul!(uxy, K₁mc, x)
    
    # x' = vx
    dx[1] = x[3]
    # y' = vy
    dx[2] = x[4]
    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    dx[3] = (n2*x[1] + two_n*x[4]) + ((μ_r2 - μ_rc3*x_1) + uxy[1])
    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    dx[4] = (n2*x[2] - two_n*x[3]) - (μ_rc3*x[2] - uxy[2])
    # t' = 1
    dx[5] = one(t)

    return dx
end

# dynamics in the 'rendezvous attempt' mode
@taylorize function spacecraft_rendezvous_attempt!(t, x, dx)
    x_1 = r + x[1]
    x_12 = x_1^2
    x_22 = x[2]^2
    rc = sqrt(x_12 + x_22)
    rc3 = rc^3
    μ_rc3 = μ / rc3

    uxy = Array{typeof(x[1])}(undef, size(K₂mc,1))
    mymul!(uxy, K₂mc, x)

    # x' = vx
    dx[1] = x[3]
    # y' = vy
    dx[2] = x[4]
    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    dx[3] = (n2*x[1] + two_n*x[4]) + ((μ_r2 - μ_rc3*x_1) + uxy[1])
    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    dx[4] = (n2*x[2] - two_n*x[3]) - (μ_rc3*x[2] - uxy[2])
    # t' = 1
    dx[5] = one(t)

    return dx
end

# dynamics in the 'aborting' mode
@taylorize function spacecraft_aborting!(t, x, dx)
    x_1 = r + x[1]
    x_12 = x_1^2
    x_22 = x[2]^2
    rc = sqrt(x_12 + x_22)
    rc3 = rc^3
    μ_rc3 = μ / rc3

    # x' = vx
    dx[1] = x[3]
    # y' = vy
    dx[2] = x[4]
    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x)
    dx[3] = (n2*x[1] + two_n*x[4]) + (μ_r2 - μ_rc3*x_1)
    # vy' = n²y - 2n*vx - μ/(rc^3)y
    dx[4] = (n^2*x[2] - 2*n*x[3]) - μ/(rc^3)*x[2]
    dx[4] = (n2*x[2] - two_n*x[3]) - μ_rc3*x[2]

    # t' = 1
    dx[5] = one(t)

    return dx
end
