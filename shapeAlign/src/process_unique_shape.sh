#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=14
#SBATCH -p largenode

# Get shape parameters for subsets of unique sites

TFS=(Abf1)
# TFS=(Abf1 Rap1 Reb1)
CLASSES=(Fast Slow)
SHAPES=(HelT MGW ProT Roll)

SHAPE_ALN=~/scripts/shapeAlign2/bin/shapeAlign

# Alignment start and end
ST=75
EN=125

# Max and min align shift
MINS=-25
MAXS=25

# Positions to ignore
IST=97
IEN=102

PARAMS="-start ${ST} -end ${EN} -min ${MINS} -max ${MAXS} -istart ${IST} -iend ${IEN}"


for C in ${CLASSES[@]}; do    # Loop over site classes
    echo "Processing: ${C}"
    for CDIR in `ls -d ${C}*`; do # Loop over each directory for site class
        for TF in ${TFS[@]}; do
            echo -e "\t${CDIR}\t${TF}"
            for S in ${SHAPES[@]}; do # Loop over each shape param
                grep -Ff ${CDIR}/${TF}*.acc ../../shape/${TF}/chec_sites/${TF}*${S}*.tab > ${CDIR}/${TF}.${S}.tab
            done

            cd ${CDIR}
            awk '{print $1}' ${TF}.HelT.tab > ${TF}.names
            ${SHAPE_ALN} ${PARAMS} -n ${TF}.names -f ${TF}*.tab
            cd ..
        done
    done
done


for S in ${SHAPES[@]}; do
    grep -Ff Rand/rand.acc /home/skasinat/scripts/chec_reanalysis/121515_max_centered/unique_shape/Rap1_rand_12_8/*${S}*.tab > Rand/rand.${S}.tab
done
cd Rand
awk '{print $1}' rand.HelT.tab > rand.names
${SHAPE_ALN} ${PARAMS} -n rand.names -f rand*.tab
cd ..
