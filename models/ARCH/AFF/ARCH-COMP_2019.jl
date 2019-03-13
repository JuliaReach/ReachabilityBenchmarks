# ISS benchmark
include("SpaceStation/iss_benchmark.jl")

# Spacecraft Rendezvous benchmark
include("Rendezvous/SpacecraftRendezvous.jl")
system, options = spacecraft()
res = run_spacecraft(system, options)
@assert res.satisfied
@time run_spacecraft(system, options)

# Powertrain benchmark
# N/A

# Building benchmark
include("Building/building_benchmark.jl")

# Platooning benchmark
include("Platooning/Platooning.jl")
system, options = platooning()
res = run_platooning(system, options)
@assert res.satisfied
@time run_platooning(system, options)

# Gearbox benchmark
include("Gearbox/Gearbox.jl")
system, options = gearbox()
res = run_gearbox(system, options)
@assert res.satisfied
@time run_gearbox(system, options)

nothing
