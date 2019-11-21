using Dates, ReachabilityBenchmarks

# print current time
now()

# ISS benchmark
include("SpaceStation/iss_benchmark.jl")

# Spacecraft Rendezvous benchmark
include("Rendezvous/SpacecraftRendezvous_benchmark.jl")

# Powertrain benchmark
# N/A

# Building benchmark
include("Building/building_benchmark.jl")

# Platooning benchmark
include("Platooning/Platooning_benchmark.jl")

# Gearbox benchmark
# include("Gearbox/Gearbox.jl")
# system, options = gearbox()
# res = run_gearbox(system, options)
# @assert res.satisfied
# @time run_gearbox(system, options)

nothing
