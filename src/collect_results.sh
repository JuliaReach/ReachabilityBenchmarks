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

function convert_to_scientific() {
    if [ -z $(echo "$1" | grep "e") ]; then
        # decimal notation
        I=$(echo $1 | awk '{print int($0)}')
        D=$(echo $1 | awk '{print $0-int($0)}' | sed 's/0\.\([0-9]*\)/\1/')
        if [ $I -eq 0 ]; then
            # number is 0.XXX
            EXPONENT=-1
            while [ ${D:0:1} -eq 0 ]; do
                D=${D:1}
                EXPONENT=$(($EXPONENT - 1))
            done
            L_D=${#D}
            if [ $L_D -eq 1 ]; then
                # one digit remains
                AFTER="${D:1:1}0"
            else
                # at least two digits remain
                AFTER="${D:1:2}"
            fi
            BEFORE="${D:0:1}"
        else
            # number is bigger than 1
            L_I=${#I}
            if [ $L_I -eq 1 ]; then
                # one-digit integer part
                AFTER="${D:0:2}"
                EXPONENT=0
            elif [ $L_I -eq 2 ]; then
                # two-digit integer part
                AFTER="${I:1:1}${D:0:1}"
                EXPONENT=1
            else
                # at least three-digit integer part
                AFTER=${I:1:2}
                EXPONENT=$(($L_I - 1))
            fi
            BEFORE="${I:0:1}"
        fi
        if [ $EXPONENT -eq 0 ]; then
            SHARED="$BEFORE.$AFTER"
        else
            SHARED="$BEFORE.${AFTER}e$EXPONENT"
        fi
    else
        # scientific notation
        E=$(echo $1 | sed 's/.*e\(.*\)/\1/')
        # take first 4 characters
        N=$(echo $1 | sed 's/\(....\).*e.*/\1/')
        SHARED="${N}e$E"
    fi
}

function collect() {
    for MODE in ${MODES[@]}; do
        OUTPUT_FILE="res-$MODE"
        
        if [ -f "$OUTPUT_FILE" ]; then
            echo "File $OUTPUT_FILE already exists, removing it..."
            rm "$OUTPUT_FILE"
        fi
        touch "$OUTPUT_FILE"
        
        echo -e "Model\tdiscretize\tanalyze\t\ttotal\n" >> "$OUTPUT_FILE"
        for MODEL in ${MODELS[@]}; do
            FILE="$MODEL-$MODE"
            
            if [ ! -f "$FILE" ]; then
                continue
            fi
            TIMES=$(grep -o "elapsed time: [0-9]*\.[0-9]*\(e-[0-9]\)\? seconds$" "$FILE" \
                  | grep -o "[0-9]*\.[0-9]*\(e-[0-9]\)\?*")
            
            # disretization
            COUNTER=1
            SHARED=""
            for T in $TIMES; do
                if [ $COUNTER -eq 1 ]; then
                    DISC=$T
                    convert_to_scientific $T
                    DISC_SC=$SHARED
                elif [ $COUNTER -eq 4 ]; then
                    # total runtime has to include discretization time
                    ANALYZE=$T
                    convert_to_scientific $ANALYZE
                    ANALYZE_SC=$SHARED
                    DISC_DEC=$(echo "$DISC" | awk '{ print +$1 }')
                    ANALYZE_DEC=$(echo "$ANALYZE" | awk '{ print +$1 }')
                    TOTAL=$(echo "$DISC_DEC + $ANALYZE_DEC" | bc | awk '{ print +$1 }')
                    convert_to_scientific $TOTAL
                    TOTAL_SC=$SHARED
                fi
                COUNTER=$(($COUNTER + 1))
            done
            
            echo -e "$MODEL\t$DISC\t$ANALYZE\t$TOTAL\n\t$DISC_SC\t\t$ANALYZE_SC\t\t$TOTAL_SC" >> "$OUTPUT_FILE"
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
#     touch $RESULT_FILE
# fi
# 
# for EX in ${NAMES[@]}; do
#     cd $EX
#     if [ ! -e $DATA_FILE ]; then
#         echo "no $DATA_FILE file found, skipping $EX..."
#         cd ..
#         continue
#     fi
#     echo ---$'\n'Results for model $EX >> "../$RESULT_FILE"
#     
#     for CATEGORY in "${CATEGORIES[@]}"; do
#         TIME=$(grep -A 1 "$CATEGORY" $DATA_FILE | tail -1 | sed -e 's/elapsed time: //' -e 's/ seconds//')
#         echo "$CATEGORY: $TIME" >> "../$RESULT_FILE"
#     done
#     echo ---$'\n' >> "../$RESULT_FILE"
#     cd ..
# done
