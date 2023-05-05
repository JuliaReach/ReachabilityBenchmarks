#!/bin/bash
#
# Example commands:
# 1) compute() for N and δ given in array above:
#   COMMAND="include(\"$EX.jl\"); compute(); compute(${STEPS[$i]}, ${DELTAS[$i]}, $OPTIONS)"
# 2) tune_δ():
#   COMMAND="include(\"$EX.jl\"); tune_δ(compute, ${T_HORIZ[$i]})"

# model names
# NAMES=(   "crane" "motor" "helicopter" "building" "pde"  "cdplayer" "heat" "iss"  "beam" "mna1" "fom" )
# # number of time steps
# STEPS=(   1500    20000   400          10000      6667   100        20000  2000   40000  20000   4000  )
# # STEPS=(   10      10      10           10         10     10         10     10     10     10     10  )
# # size of time step; -1. = use default step size in model file
# DELTAS=(  -1.     -1.     -1.          -1.        -1.    -1.        -1.    -1.    -1.    -1.    -1. )
# # time horizon
# T_HORIZ=( 15.     20.     20.          20.        20.    2.         20.    20.    20.    20.    20. )
# # options
# OPTIONS="spaceExOptions" # "hylaaOptions"

# fixed delta, SLICOT only without MNA5
# NAMES=( "motor" "building" "pde"  "heat" "iss"  "beam" "mna1" "fom" )
# STEPS=(   20000   20000      20000  20000  20000  20000  20000  20000  )
# DELTAS=(  1e-3    1e-3       1e-3   1e-3   1e-3   1e-3   1e-3   1e-3 )
# T_HORIZ=( 20.     20.        20.    20.    20.    20.    20.    20. )
# OPTIONS="spaceExOptions"

# property checking, SLICOT only without PDE/ISS/FOM
NAMES=(   "motor" "building" "heat" "beam" "mna1" "mna5" )
STEPS=(   4000    4000       4000   4000   4000   4000   )
DELTAS=(  0.005   0.005      0.005  0.005  0.005  0.005  )
OPTIONS="hylaaOptions"

echo "Running ${#NAMES[@]} benchmarks..."
i=0
for EX in ${NAMES[@]}; do
	echo -----------------------
	cd $EX
	COMMAND="include(\"$EX.jl\"); compute(3, ${DELTAS[$i]}, $OPTIONS); compute(${STEPS[$i]}, ${DELTAS[$i]}, $OPTIONS)"
# 	COMMAND="include(\"$EX.jl\"); tune_δ(compute, 20.)"
	echo $(($i+1))": Running Julia on $EX example with the following command"
	echo "   julia -O3 -e '$COMMAND'"
 	julia -O3 -e "$COMMAND" > out.txt # remove '> out.txt' for seeing the output
	cd ..
	i=$i+1
done

# for EX in ${NAMES[@]}; do
# 	echo -----------------------
# 	cd $EX
# 	echo $(($i+1))": Running SpaceEx (2D) on $EX example"
# 	./spaceex.sh "--sampling-time 0.001"
# 	cd ..
# 	i=$i+1
# done

# for EX in ${NAMES[@]}; do
# 	echo -----------------------
# 	cd $EX
# 	echo $(($i+1))": Running SpaceEx (full) on $EX example"
# 	./spaceex.sh "--sampling-time 0.001 --directions box"
# 	cd ..
# 	i=$i+1
# done

# OPTIONSII="spaceExOptions"
# echo "Running ${#NAMES[@]} benchmarks..."
# i=0
# for EX in ${NAMES[@]}; do
# 	echo -----------------------
# 	cd $EX
# 	COMMAND="include(\"$EX.jl\"); compute(3, ${DELTAS[$i]}, $OPTIONSII); compute(${STEPS[$i]}, ${DELTAS[$i]}, $OPTIONSII)"
# # 	COMMAND="include(\"$EX.jl\"); tune_δ(compute, ${T_HORIZ[$i]})"
# 	echo $(($i+1))": Running Julia on $EX example with the following command"
# 	echo "   julia -O3 -e '$COMMAND'"
#  	julia -O3 -e "$COMMAND" > out.txt # remove '> out.txt' for seeing the output
# 	cd ..
# 	i=$i+1
# done
