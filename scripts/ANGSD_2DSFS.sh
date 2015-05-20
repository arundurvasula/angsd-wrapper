#!/usr/bin/env bash

set -e
set -u

#   Source the common configuration file
source scripts/common.conf
# load utils functions
source ${SCRIPTS_DIR}/utils.sh

#Defaults
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
OVERRIDE=false

# load variables from supplied config file
load_config $1

# calculated variables
N_IND1=`wc -l < ${TAXON_LIST1}`
N_CHROM1=`expr 2 \* ${N_IND1}`
N_IND2=`wc -l < ${TAXON_LIST2}`
N_CHROM2=`expr 2 \* ${N_IND2}`

# For 1st taxon
if file_exists "${RESULTS_DIR}/${TAXON1}_Intergenic.mafs.gz" && [ "$OVERRIDE" = "false" ]; then 
    echo "maf already exists and OVERRIDE=false, skipping angsd -bam..."
elif [[ ${REGIONS} == */* ]]; then
    ${ANGSD_DIR}/angsd\
        -bam ${TAXON_LIST1}\
        -out ${RESULTS_DIR}/${TAXON1}_Intergenic\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING1}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND1}\
        -minInd ${MIN_IND1}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -rf ${REGIONS}
else
    ${ANGSD_DIR}/angsd\
        -bam ${TAXON_LIST1}\
        -out ${RESULTS_DIR}/${TAXON1}_Intergenic\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING1}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND1}\
        -minInd ${MIN_IND1}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -r ${REGIONS}
fi

${ANGSD_DIR}/misc/realSFS\
    ${RESULTS_DIR}/${TAXON1}_Intergenic.saf\
    ${N_CHROM1}\
    -P ${N_CORES}\
    > ${RESULTS_DIR}/${TAXON1}_Intergenic.sfs 

# For 2nd taxon:
if [ -e "${RESULTS_DIR}/${TAXON2}_Intergenic.mafs.gz" ] && [ "$OVERRIDE" = "false" ]; then
    echo "maf already exists and OVERRIDE=false, skipping angsd -bam..."
elif [[ ${REGIONS} == */* ]]; then
    ${ANGSD_DIR}/angsd\
        -bam ${TAXON_LIST2}\
        -out ${RESULTS_DIR}/${TAXON2}_Intergenic\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING2}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND2}\
        -minInd ${MIN_IND2}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -rf ${REGIONS}
else
    ${ANGSD_DIR}/angsd\
        -bam ${TAXON_LIST2}\
        -out ${RESULTS_DIR}/${TAXON2}_Intergenic\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING1}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND2}\
        -minInd ${MIN_IND2}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -r ${REGIONS}
fi

${ANGSD_DIR}/misc/realSFS\
    ${RESULTS_DIR}/${TAXON2}_Intergenic.saf\
    ${N_CHROM2}\
    -P ${N_CORES}\
    > ${RESULTS_DIR}/${TAXON2}_Intergenic.sfs

#find intersecting regions
gunzip -c ${RESULTS_DIR}/${TAXON1}_Intergenic.saf.pos.gz ${RESULTS_DIR}/${TAXON2}_Intergenic.saf.pos.gz | sort -S 50% | uniq -d | sort -k1,1 -S 50% > ${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt

#index intersect file
${ANGSD_DIR}/angsd sites index ${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt

# calculate allele frequencies only on sites in both populations
if [[ ${REGIONS} == */* ]]; then
    ${ANGSD_DIR}/angsd\
        -bam ${TAXON_LIST1}\
        -out ${RESULTS_DIR}/${TAXON1}_Intergenic_Conditioned\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING1}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND1}\
        -minInd ${MIN_IND1}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -rf ${REGIONS}\
        -sites ${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt
else
    ${ANGSD_DIR}/angsd\
        -bam ${TAXON_LIST1}\
        -out ${RESULTS_DIR}/${TAXON1}_Intergenic_Conditioned\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING1}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND1}\
        -minInd ${MIN_IND1}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -r ${REGIONS}\
        -sites ${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt
fi

if [[ ${REGIONS} == */* ]]; then
    ${ANGSD_DIR}/angsd\
        -bam ${TAXON_LIST2}\
        -out ${RESULTS_DIR}/${TAXON2}_Intergenic_Conditioned\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING2}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND2}\
        -minInd ${MIN_IND2}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -rf ${REGIONS}\
        -sites ${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt
else
    ${ANGSD_DIR}/angsd\
        -bam ${TAXON_LIST2}\
        -out ${RESULTS_DIR}/${TAXON2}_Intergenic_Conditioned\
        -doMajorMinor ${DO_MAJORMINOR}\
        -doMaf ${DO_MAF}\
        -indF ${TAXON_INBREEDING2}\
        -doSaf ${DO_SAF}\
        -uniqueOnly ${UNIQUE_ONLY}\
        -anc ${ANC_SEQ}\
        -minMapQ ${MIN_MAPQ}\
        -minQ ${MIN_BASEQUAL}\
        -nInd ${N_IND2}\
        -minInd ${MIN_IND2}\
        -baq ${BAQ}\
        -ref ${REF_SEQ}\
        -GL ${GT_LIKELIHOOD}\
        -P ${N_CORES}\
        -rf ${REGIONS}\
        -sites ${RESULTS_DIR}/intersect.${TAXON1}.${TAXON2}_intergenic.txt
fi

# estimate joint SFS using realSFS
${ANGSD_DIR}/misc/realSFS 2dsfs\
    ${RESULTS_DIR}/${TAXON1}_Intergenic_Conditioned.saf ${RESULTS_DIR}/${TAXON2}_Intergenic_Conditioned.saf\
    ${N_IND1}\
    ${N_IND2}\
    -P ${N_CORES}\
    > ${RESULTS_DIR}/2DSFS_Intergenic.${TAXON1}.${TAXON2}.sfs
