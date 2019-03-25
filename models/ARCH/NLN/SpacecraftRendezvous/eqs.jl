const μ = 3.986e14 * 60^2
const r = 42164e3
const mc = 500.0
const n = sqrt(μ / r^3)

#K₁ = [-28.8287 0.1005 -1449.9754 0.0046;
#      -0.087 -33.2562 0.00462 -1451.5013]
const K₁11 = -28.8287; const K₁12 = 0.1005; const K₁13 = -1449.9754; const K₁14 = 0.0046;
const K₁21 = -0.087; const K₁22 = -33.2562; const K₁23 = 0.00462; const K₁24 = -1451.5013;

#K₂ = [-288.0288 0.1312 -9614.9898 0.0;
#      -0.1312 -288.0 0.0 -9614.9883]
const K₂11 = -288.0288; const K₂12 = 0.1312; const K₂13 = -9614.9898; const K₂14 = 0.0;
const K₂21 = -0.1312; const K₂22 = -288.0; const K₂23 = 0.0; const K₂24 = -9614.9883;

# dynamics in the 'approaching' mode
@taylorize function spacecraft_approaching!(t, x, dx)

    rc = sqrt((r + x[1])^2 + x[2]^2)
    uxy1 = (K₁11 * x[1] + K₁12 * x[2]) + (K₁13 * x[3] + K₁14 * x[4])
    uxy2 = (K₁21 * x[1] + K₁22 * x[2]) + (K₁23 * x[3] + K₁24 * x[4])

    # x' = vx
    dx[1] = x[3]

    # y' = vy
    dx[2] = x[4]

    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    dx[3] = (n^2*x[1] + 2*(n*x[4])) + ((μ/(r^2) - μ/(rc^3)*(r + x[1])) + uxy1/mc)

    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    dx[4] = (n^2*x[2] - 2*(n*x[3])) - (μ/(rc^3)*x[2] - uxy2/mc)

    # t' = 1
    dx[5] = 1.0 + 0.0*x[1]

    return dx
end

# dynamics in the 'rendezvous attempt' mode
@taylorize function spacecraft_rendezvous_attempt!(t, x, dx)

    rc = sqrt((r + x[1])^2 + x[2]^2)
    uxy1 = (K₂11 * x[1] + K₂12 * x[2]) + (K₂13 * x[3] + K₂14 * x[4])
    uxy2 = (K₂21 * x[1] + K₂22 * x[2]) + (K₂23 * x[3] + K₂24 * x[4])

    # x' = vx
    dx[1] = x[3]

    # y' = vy
    dx[2] = x[4]

    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    dx[3] = (n^2*x[1] + 2*(n*x[4])) + ((μ/(r^2) - μ/(rc^3)*(r + x[1])) + uxy1/mc)

    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    dx[4] = (n^2*x[2] - 2*(n*x[3])) - (μ/(rc^3)*x[2] - uxy2/mc)

    # t' = 1
    dx[5] = 1.0 + 0.0*x[1]

    return dx
end

# dynamics in the 'aborting' mode
@taylorize function spacecraft_aborting!(t, x, dx)

    rc = sqrt((r + x[1])^2 + x[2]^2)

    # x' = vx
    dx[1] = x[3]

    # y' = vy
    dx[2] = x[4]

    # vx' = n²x + 2n*vy + μ/(r^2) + μ/(rc^3)*(r+x)
    dx[3] = (n^2*x[1] + 2*n*x[4]) + μ/(r^2) + μ/(rc^3)*(r + x[1])

    # vy' = n²y - 2n*vx - μ/(rc^3)y
    dx[4] = (n^2*x[2] - 2*n*x[3]) - μ/(rc^3)*x[2]

    # t' = 1
    dx[5] = 1.0 + 0.0*x[1]

    return dx
end
