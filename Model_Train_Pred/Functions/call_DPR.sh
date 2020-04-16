#!/bin/bash

file_dir=$1
dpr_num=$2
TargetID=$3

DPR=`pwd`/Model_Train_Pred/Functions/DPR

cd ${file_dir}

BIMBAM=`pwd`/bimbam
PHENO=`pwd`/pheno
SNP=`pwd`/SNP_annot

echo $TargetID

${DPR} \
-g ${BIMBAM}/${TargetID}_bimbam.txt \
-p ${PHENO}/${TargetID}_pheno.txt \
-a ${SNP}/${TargetID}_snp_annot.txt \
-dpr ${dpr_num} \
-o DPR_${TargetID}


