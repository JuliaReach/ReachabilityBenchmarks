r"""

Transform a linear combination into a matrix in txt or MAT format.

This script can be used to transform input typically from the SX format
(SpaceEx models), to a projection matrix that is read in Julia.

USAGE:

1. Copy the linear combination as a string, and store it in the
format: ``'y = 4.3036e-9*x136 - 0.000044272*x137 - 2.0608e-14*x138'``
The filename extension should be ``'.sage'``.

2. Call this Sage script as follows:

        $ sage outProjectionMatrix.sage n input.sage output_format

    where:

     - ``n``             -- is an integer with the dimension of the system
     - ``input.sage``    -- is the file which conatins the linear combination
     - ``output_format`` -- can be either ``txt`` or ``MAT``

3. The result is recovered in the file ``'out'``, whose filename extension is
given by the choice of ``output_format``.

EXAMPLES:

    We create a linear combination of a system in 270 variables::

        $ echo y =  -6.314e-4*x1 + 4.3036e-3*x136 - 0.000044272*x137 - 2.0608e-14*x138 >> y.sage
        $ cat y.sage
        y = -6.314e-4*x1 + 4.3036e-3*x136 - 0.000044272*x137 - 2.0608e-14*x138
        $ sage outProjectionMatrix.sage 270 y.sage txt

    We have produced a projection matrix ``M``::

        $ cat out.txt
        M[1, 1] = -0.000631400000000000
        M[1, 136] = 0.00430360000000000
        M[1, 137] = -0.0000442720000000000
        M[1, 138] = -2.06080000000000e-14

    Assuming that the file ``pde.sage`` contains a linear combination (see below),
    and that the model is 84-dimensional::

        $ sage outProjectionMatrix.sage 84 pde.sage MAT
        $ julia
        julia> using MAT; M = sparse(read(matopen("out.mat"), "M"))
        1x84 SparseMatrixCSC{Float64,Int64} with 83 stored entries:
          [1 ,  1]  =  9.3391
          [1 ,  2]  =  2.3546
          [1 ,  3]  =  7.5703
          [1 ,  4]  =  2.3216
          ...
          [1 , 80]  =  7.6368
          [1 , 81]  =  0.66947
          [1 , 82]  =  2.4018
          [1 , 83]  =  3.5422
          [1 , 84]  =  6.0766

    Input for the second example (``'pde.sage'``):

    y = 7.6248*x1 + 2.3546*x2 + 7.5703*x3 + 2.3216*x4 + 4.6086*x5 + 0.3538*x6 + 5.2422*x7 + 3.9923*x8 + 8.9943*x9 + 1.7143*x1 + 0.24718*x11 + 5.9218*x12 + 5.7022*x13 + 5.8246*x14 + 3.4134*x15 + 5.5956*x16 + 3.9767*x17 + 9.4569*x18 + 8.765*x19 + 9.8845*x20 + 3.9137*x21 + 9.2012*x22 + 4.5751*x23 + 6.6164*x24 + 5.2145*x25 + 9.829*x26 + 8.9967*x27 + 6.302*x28 + 3.4683*x29 + 3.0584*x30 + 9.0124*x31 + 1.066*x32 + 2.9632*x33 + 8.5156*x34 + 7.2941*x35 + 6.4278*x36 + 2.6819*x37 + 7.2813*x38 + 8.922*x39 + 5.0993*x40 + 5.5309*x41 + 8.8087*x42 + 8.5626*x43 + 5.9954*x44 + 4.8225*x45 + 0.12009*x46 + 2.7879*x47 + 5.7413*x48 + 8.2792*x49 + 7.0457*x50 + 3.4186*x51 + 4.0203*x52 + 5.389*x53 + 1.8974*x54 + 4.3469*x55 + 4.0868*x56 + 2.0564*x57 + 2.9339*x58 + 2.5612*x59 + 3.1261*x60 + 5.4868*x61 + 8.3082*x62 + 2.975*x63 + 8.7716*x64 + 7.5489*x65 + 3.2524*x66 + 2.5812*x67 + 7.3413*x68 + 4.3744*x69 + 3.1246*x70 + 9.0366*x71 + 6.9508*x72 + 7.5679*x73 + 2.6664*x74 + 6.1122*x75 + 1.7228*x76 + 1.5188*x77 + 6.8394*x78 + 2.1795*x79 + 7.6368*x80 + 0.66947*x81 + 2.4018*x82 + 3.5422*x83 + 6.0766*x84

NOTES:

Use the `{1:.8}` parameter to truncate the number of decimal digits in the coefficients.
Without the double dot modifier, it prints the full number (recommended).

AUTHORS:

- Marcelo Forets (2017-09-06)
"""
import sys, os
if len(sys.argv) < 4:
    raise Exception("wrong number of input arguments. usage: $ sage outProjectionMatrix.sage n input.sage txt|MAT")

n = int(sys.argv[1])
var(['x'+str(1+i) for i in range(n)]);  # inject x1, ..., xn into global namespace
load(sys.argv[2])
output_format = sys.argv[3]

if output_format == "txt":
    print("writing projection matrix M in out.txt")
    with open("out.txt", 'w') as fp:
        for v in y.variables():
            print >> fp, "M[1, {0}] = {1}".format(int(str(v)[1:]), y.coefficient(v));

elif output_format == "MAT":
    print("writing projection matrix M in out.MAT")
    from scipy.io import savemat
    M = zero_matrix(RDF, 1, n, sparse=True)
    for v in y.variables():
        M[0, int(str(v)[1:])-1] = y.coefficient(v)
    savemat("out", {"M": M})

# remove auxiliary .sage.py file (preparsing temp)
os.remove("outProjectionMatrix.sage.py")
