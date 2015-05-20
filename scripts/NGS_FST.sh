#!/usr/bin/env bash

set -e
set -u

#   source common config file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

#Need to have run ANGSD_2DSFS or have files in correct name format

DO_SAF=2
UNIQUE_ONLY=0
MIN_BASEQUAL=20
BAQ=1
MIN_IND1=1
MIN_IND2=1
GT_LIKELIHOOD=2
MIN_MAPQ=30
N_CORES=16
DO_MAJORMINOR=1
DO_MAF=1
BLOCK_SIZE=20000

load_config $1

TAXON1_LIST=${DATA_DIR}/${TAXON1}_samples.txt
TAXON2_LIST=${DATA_DIR}/${TAXON2}_samples.txt
POP1_SFS=${RESULTS_DIR}/${TAXON1}_Intergenic_Conditioned
POP2_SFS=${RESULTS_DIR}/${TAXON2}_Intergenic_Conditioned
INTERSECT=${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt
2DSFS=${RESULTS_DIR}/2DSFS_Intergenic.${TAXON1}.${TAXON2}.sfs

N_IND_1=`wc -l < ${TAXON1_LIST}`
N_IND_2=`wc -l < ${TAXON2_LIST}`

#check for sfs, intersect file, and 2dsfs from ANGSD_2DSFS
#exit with error if any don't exist
if [ -e "${POP1_SFS}"]; then 
    echo "saf for ${TAXON1} exists, continuing to check for ${TAXON2} saf..."
else echo "saf for ${TAXON1} does not exist, exiting..." >&2; exit 1
fi

if [ -e "${POP2_SFS}"]; then 
    echo "saf for ${TAXON2} exists, continuing to check for intersect..."
else echo "saf for ${TAXON2} does not exist, exiting..." >&2; exit 1 
fi

if [ -e "${INTERSECT}"]; then 
    echo "intersect exists, continuing to check for 2dsfs..."
else echo "intersect does not exist, exiting..." >&2; exit 1
fi

if [ -e "${2DSFS}"]; then 
    echo "2dsfs exists, continuing to analysis..."
else echo "2dsfs does not exist, exiting..." >&2; exit 1
fi

#make sure intersect index exists
#if not, make it
if [ ! -e "${INTERSECT}.idx"]; then
    ${ANGSD_DIR}/angsd sites index ${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt
fi

# get number of sites and individuals
N_SITES=`wc -l ${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt | cut -f 1 -d " "`

# convert ANGSD 2DSFS for ngsPopGen use
${ANGSD_DIR}/scripts/convertSFS.R ${RESULTS_DIR}/2DSFS_Intergenic.${TAXON1}.${TAXON2}.sfs\
    > ${RESULTS_DIR}/converted.2DSFS_Intergenic.${TAXON1}.${TAXON2}.sfs

# get FST
${NGS_POPGEN_DIR}/ngsFST\
    -postfiles ${POP1_SFS} ${POP2_SFS}\
    -priorfile ${RESULTS_DIR}/converted.2DSFS_Intergenic.${TAXON1}.${TAXON2}.sfs\
    -nind ${N_IND_1} ${N_IND_2}\
    -nsites ${N_SITES}\
    -block_size ${BLOCK_SIZE}\
    -n_threads ${N_CORES}\
    -outfile ${RESULTS_DIR}/${TAXON1}.${TAXON2}.fst
