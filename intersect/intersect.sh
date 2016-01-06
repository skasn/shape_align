#!/bin/bash

# Given a set of BED files containing single base coordinates,
# perform intersections at a variety of window widths

# Start and end widths: these are half window widths
# i.e., LST=1 corresponds to 1 bp upstream and 1 bp downstream
LST=1
LEN=500

# Window length step size
LSTEP=1

# Chrom sizes
SC3_SIZES=~/genomes/sacCer3.sizes

# Function for returning comma-separated or tab-delimited
# array values
function csvout { local IFS=","; echo "$*"; }
function tabout { local IFS=$'\t'; echo "$*"; }

module load bedtools

for ((w=${LST}; w<=${LEN}; w+=${LSTEP})); do
	W2=$((2*w))
	echo -e "Win size: $W2"
	
	CTR=0
	SLOPPED=()	# Names of windowed files
	for arg; do
		SUFF=win.${W2}.${arg##*.}
		SLOPPED[${CTR}]=${arg%.*}.${SUFF}
		bedtools slop -b ${W2} -i ${arg} -g ${SC3_SIZES} > ${SLOPPED[${CTR}]}
		CTR=$((CTR+1))
	done	
	
	NUNIQ=()	# Array for counting unique entries in each file
	# Loop over all input files
	for  ((i=1; i<=$#; i++)); do
		# Make a copy of the windowed file names array
		MYARGS=("${SLOPPED[@]}")
		A=${MYARGS[$((i-1))]}
		
		# Remove A from MYARGS
		MYARGS=("${MYARGS[@]::$((i-1))}" "${MYARGS[@]:${i}}")
		B="${MYARGS[@]}"
		
		# Perform intersection
		NUNIQ[$i]=`bedtools intersect -v -a ${A} -b ${B} | wc -l`
	done

	tabout ${w} ${NUNIQ[@]} >> isect.txt
	
	# Clean up
	for SLOPFILE in ${SLOPPED[@]}; do
		rm ${SLOPFILE}
	done

done

