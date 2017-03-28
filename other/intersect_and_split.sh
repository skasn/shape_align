#!/bin/bash

# Given a set of BED files containing single base coordinates
# and a window width, get unique sites for each ChEC class

# Pass a list of BED files

# Chrom sizes
SC3_SIZES=sacCer3.sizes

# Classifications of sites (must be unique and in BED file)
CLASS=(Fast Slow)

# Random sites file
RSITES=random.bed

# Function for returning comma-separated or tab-delimited
# array values
function csvout { local IFS=","; echo "$*"; }
function tabout { local IFS=$'\t'; echo "$*"; }

# Window half-widths
HWIDTHS=(25 50 75 100 125 150 175 200 225 250 275 300)

for W in ${HWIDTHS[@]}; do
    W2=$((2*W))
    echo -e "Win size: $W2"

    CTR=0
    SLOPPED=()  # Names of windowed files
    for arg; do
        SUFF=win.${W2}.${arg##*.}
        SLOPPED[${CTR}]=`basename ${arg%.*}.${SUFF}`
        bedtools slop -b ${W} -i ${arg} -g ${SC3_SIZES} > ${SLOPPED[${CTR}]}
        CTR=$((CTR+1))
    done

    NUNIQ=()    # Array for counting unique entries in each file
    # Loop over all input files
    for  ((i=1; i<=$#; i++)); do
        # Make a copy of the windowed file names array
        MYARGS=("${SLOPPED[@]}")
        A=${MYARGS[$((i-1))]}

        # Remove A from MYARGS
        MYARGS=("${MYARGS[@]::$((i-1))}" "${MYARGS[@]:${i}}")
        B="${MYARGS[@]}"

        # Perform intersection
        bedtools intersect -v -a ${A} -b ${B} > ${A%.*}.unique.bed
        NUNIQ[$i]=`wc -l ${A%.*}.unique.bed | awk '{print $1}'`

        for C in ${CLASS[@]}; do
            grep ${C} ${A%.*}.unique.bed | awk '{print $4}' > ${A%.*}.unique.${C}.w${W2}.acc
        done
        rm ${A%.*}.unique.bed
    done

    # Print the intersection info
    tabout ${W} ${NUNIQ[@]} >> isect.win.${W2}.txt

    # Get the unique random sites
    RSLOP=${RSITES%.*}.win${W2}.bed
    RUNIQ=${RSITES%.*}.w${W2}.acc
    bedtools slop -b ${W2} -i ${RSITES} -g ${SC3_SIZES} > ${RSLOP}

    bedtools intersect -v -a ${RSLOP} -b ${SLOPPED[@]} > ${RUNIQ}


    # Clean up
    for SLOPFILE in ${SLOPPED[@]}; do
        rm ${SLOPFILE}
    done
    rm ${RSLOP}
done
