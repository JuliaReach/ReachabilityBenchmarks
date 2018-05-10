#!/bin/bash

## Settings 
# model names
MODELS=("motor" "building" "pde" "heat" "iss" "beam" "mna1" "fom" "mna5")
# individual file name suffix
MODES=( \
       "reach-1D-fixedstep-onevar.txt" \
       "reach-1D-fixedstep-allvars.txt" \
       "check-1D-dense-givenstep.txt" \
       "check-1D-discrete-fixedstep.txt" \
       "reach-2D-box-fixedstep-twovars.txt" \
       "reach-2D-eps-fixedstep-twovars.txt" \
      )

function collect(){
    for MODE in ${MODES[@]}; do
        OUTPUT_FILE="res-$MODE"
        
        if [ -f "$OUTPUT_FILE" ]; then
            echo "File $OUTPUT_FILE already exists, removing it..."
            rm "$OUTPUT_FILE"
        fi
        touch "$OUTPUT_FILE"
        
        echo -e "Model\tdiscretize\ttotal\n" >> "$OUTPUT_FILE"
        for MODEL in ${MODELS[@]}; do
            FILE="$MODEL-$MODE"
            
            if [ ! -f "$FILE" ]; then
                continue
            fi
            TIMES=$(grep -o "elapsed time: [0-9]*\.[0-9]*\(e-[0-9]\)\? seconds$" "$FILE" \
                  | grep -o "[0-9]*\.[0-9]*\(e-[0-9]\)\?*")
            
            # disretization
            COUNTER=1
            for T in $TIMES; do
                if [ $COUNTER -eq 1 ]; then
                    DISC=$T
                elif [ $COUNTER -eq 4 ]; then
                    TOTAL=$T
                fi
                COUNTER=$(($COUNTER + 1))
            done
            
            echo -e "$MODEL\t$DISC\t$TOTAL" >> "$OUTPUT_FILE"
        done
    done
}

pushd . > /dev/null
collect
echo ""
popd > /dev/null



## old script
# 
# # collects all benchmark results from the following list
# NAMES=( "crane" "motor" "helicopter" "building" "pde"  "cdplayer" "heat" "iss"  "beam" "mna1" "fom" "mna5" )
# 
# # collect all timings for the following categories
# CATEGORIES=( "System Construction" "Time Discretization" "Reachable States Computation" "Projection" )
# 
# # name of the data file
# DATA_FILE="out.txt"
# 
# # name of the result file
# RESULT_FILE="results-$(date +%Y%m%d%H%M%S).txt"
# if [ ! -e $RESULT_FILE ]; then
# 	touch $RESULT_FILE
# fi
# 
# for EX in ${NAMES[@]}; do
# 	cd $EX
# 	if [ ! -e $DATA_FILE ]; then
# 		echo "no $DATA_FILE file found, skipping $EX..."
# 		cd ..
# 		continue
# 	fi
# 	echo ---$'\n'Results for model $EX >> "../$RESULT_FILE"
# 	
# 	for CATEGORY in "${CATEGORIES[@]}"; do
# 		TIME=$(grep -A 1 "$CATEGORY" $DATA_FILE | tail -1 | sed -e 's/elapsed time: //' -e 's/ seconds//')
# 		echo "$CATEGORY: $TIME" >> "../$RESULT_FILE"
# 	done
# 	echo ---$'\n' >> "../$RESULT_FILE"
# 	cd ..
# done
