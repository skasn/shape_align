#!/bin/bash

# Overview of script:
# Input: a list of peaks (in BED format), a directory containing raw data
#        in sorted BEDGRAPH format, and window width
# Output: the same list of peaks with a temporally ordered columns with
#         the signal:noise ratio statistic
#
# Approach in outline:
# 1. Perform windowing of all peaks
# 2. Compute signal:noise ratio for each peak at each time point
# 3. Generate a table where each row in the table is a peak and columns
#    correspond to the signal:noise statistic at each time point
#
# A note on calculation of signal:noise ratio: I'm defining signal:noise
# as simple the ratio between the value at the peak maximum and the area
# under the curve (AUC) of the peak itself

# Input peak list:
BED=$1

# Directory containing raw data in sorted (i.e., "sort -k1,1 -k2,2n")
# BEDGRAPH format.
# A note on ordering of files: Pad times with zeros such that they are
# displayed in the proper temporal ordering when 'ls' is used to dislay
# the names. This is because the ordering of columns of the final statistic
# of interest is dependent on the 'ls' command based ordering.
ENDS_DIR=$2

WIN=25  # Half window width for calculating signal:noise

NULL=NA  # How to represent missing values? Represent as NA.

# For FHCRC cluster:
# module load bedtools
# SC3_SIZES=~/genomes/sacCer3.sizes   # sacCer3 genome info
# SCRIPTS=/home/skasinat/scripts      # scripts directory (bedmid.pl)

# For running on desktop
BTOOLS_PATH=~/Downloads/bedtools2-master/bin/
SCRIPTS=/Users/sivakasinathan/scripts/for_rhino/ # scripts directory
export PATH=${PATH}:${BTOOLS_PATH}:${SCRIPTS}
SC3_SIZES=/Users/sivakasinathan/Dropbox/NB/asat_inversion/chec_shape_register/sacCer3.sizes

PREFIX=`basename ${BED%.*}`
STATSBED=${PREFIX}.stats.bed  # Final output file

# Given that outputs are appended to files,
# make sure files with same names do not exist from
# prev. analyses -- if they do, get rid of them
if [ -e ${STATSBED} ]; then
  rm ${STATSBED}
fi

# Generate a file containing a single-column version of the input BED file
awk '{print $1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6 "_" $7 "_" $8 "_" $9 "_" $10 "_" $11}' ${BED} > ${PREFIX}.collapsed

echo "Getting max:AUC ratio in peaks across timepoints"

# Get midpoint, define desired window, and clean up
MID=${PREFIX}.mid.bed
bedmid.pl ${BED}
bedtools slop -b ${WIN} -i ${MID} -g ${SC3_SIZES} > ${PREFIX}.win$((WIN*2+1)).bed
rm ${MID}
BED=${PREFIX}.win$((WIN*2+1)).bed  # Redefine the BED file as the windowed BED

# Loop over all of the raw data files and compute the signal:noise
for RAWBED in `ls ${ENDS_DIR}`
do
  echo -e "\tTimepoint: ${RAWBED}"

  # Get sums across peak windows
  bedtools map -null "${NULL}" -a ${BED} -b ${ENDS_DIR}/${RAWBED} -c 4 -o sum | \
    awk '{print $12}' >  ${PREFIX}.${RAWBED}.sum

  # Get max value in peak windows
  # bedtools map -a ${BED} -b ${ENDS_DIR}/${RAWBED} -c 4 -o max | \
  #   awk '{print $12}' >  ${PREFIX}.${RAWBED}.max

  # Link the sum and max and compute the signal:noise ratio
  # paste ${PREFIX}.${RAWBED}.max ${PREFIX}.${RAWBED}.sum | \
  #   awk '{if ($2 > 0) {print $1/$2} else {print "NA"}}' \
  #   > ${PREFIX}.${RAWBED}.stats

  paste ${PREFIX}.${RAWBED}.sum > ${PREFIX}.${RAWBED}.stats

  # Clean up
  rm ${PREFIX}.${RAWBED}.max
  rm ${PREFIX}.${RAWBED}.sum
done

STATNAMES=`ls ${PREFIX}.*.stats`

# Create the final ouptut and transpose it

paste ${PREFIX}.collapsed ${STATNAMES[@]} | \
awk '{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' > ${STATSBED}

# Final clean up
rm ${PREFIX}.*.stats    # Stats files
rm ${BED}               # Windowed BED
rm ${PREFIX}.collapsed  # Single column representation of BED
# EOF
