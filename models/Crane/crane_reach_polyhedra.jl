using Polyhedra
using CDDLib
using MAT

include("../../src/Reachability/ReachabilityPolyhedra.jl")
using ReachabilityPolyhedra

function compute(algorithm="reach_homog")

# system's matrix
A = [0. 1. 0. 0. 0. 0. ;
-0.417533 -3.1931759963 39.24 0. -14.825331 11.123344 ;
0. 0. 0. 1. 0. 0. ;
0.0417533 0.31931759963 -4.905 0. 1.4825331 -1.1123344 ;
0.0638407957 -0.32473339016573 0. 0. -3.7332068901 -0.7007592976 ;
0.0853437452 -0.72366802635628 0. 0. -5.9714023436 -2.2736115136]

# number of variables (dimension)
n = size(A, 1)

# initial set of states
file = matopen("X0.mat")

#=
F = read(file, "F");
g = read(file, "g");
g = g[:];  # cast into a column vector
# this is a half-space representation. use the methods X0.A and X0.b to access
# the matrix/vectors such that X0 = {x : Ax <= b}
X0HRep = SimpleHRepresentation(F, g);
X0 = polyhedron(X0HRep, CDDLibrary());
=#

# initial states : hyperrectangle
g1 = read(file, "g1")[:]
g2 = read(file, "g1")[:]
g3 = read(file, "g1")[:]
F1 = read(file, "F1")
F2 = read(file, "F2")
F3 = read(file, "F3")

X0hat = Array{CDDPolyhedron}(3)
X0hat[1] = polyhedron(SimpleHRepresentation(F1, g1))
X0hat[2] = polyhedron(SimpleHRepresentation(F2, g2))
X0hat[3] =  polyhedron(SimpleHRepresentation(F3, g3))

# number of samples
N = 250

# ==== case with inputs === 

# this polyhedron is "null" in n dims
#nullPoly_6d = polyhedron(SimpleVRepresentation([0.0 0.0 0.0 0.0 0.0 0.0]));
# we specify the backend using a second argument, as CDDLibrary()
# otherwise, the object will be of type Polyhedra.SimplePolyhedron{6,Float64},
# and we cannot do operations with it
nullPoly_6d = polyhedron(SimpleVRepresentation(zeros(1, 6)), CDDLibrary());
U = [nullPoly_6d for i in 1:N]

# === solve ===

if algorithm == "reach_homog"
    RSet = reach_homog(A, X0hat, N)
elseif algorithm == "reach"
    RSet = reach(A, X0hat, N, U)
else
    error("algorithm unkwnown")
end



end