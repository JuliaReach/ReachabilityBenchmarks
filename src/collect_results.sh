#!/bin/bash

# collects all benchmark results from the following list
NAMES=( "crane" "motor" "helicopter" "building" "pde"  "cdplayer" "heat" "iss"  "beam" "mna1" "fom" "mna5" )

# collect all timings for the following categories
CATEGORIES=( "System Construction" "Time Discretization" "Reachable States Computation" "Projection" )

# name of the data file
DATA_FILE="out.txt"

# name of the result file
RESULT_FILE="results-$(date +%Y%m%d%H%M%S).txt"
if [ ! -e $RESULT_FILE ]; then
	touch $RESULT_FILE
fi

for EX in ${NAMES[@]}; do
	cd $EX
	if [ ! -e $DATA_FILE ]; then
		echo "no $DATA_FILE file found, skipping $EX..."
		cd ..
		continue
	fi
	echo ---$'\n'Results for model $EX >> "../$RESULT_FILE"
	
	for CATEGORY in "${CATEGORIES[@]}"; do
		TIME=$(grep -A 1 "$CATEGORY" $DATA_FILE | tail -1 | sed -e 's/elapsed time: //' -e 's/ seconds//')
		echo "$CATEGORY: $TIME" >> "../$RESULT_FILE"
	done
	echo ---$'\n' >> "../$RESULT_FILE"
	cd ..
done