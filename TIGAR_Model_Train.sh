#!/bin/bash

#####################################################################################################
# software requirement
# python 3
# tabix
###

# Variables Needed for Training
###
# model: elastic_net or DPR
# Gene_Exp: Gene Annotation and Expression Level File 
# sampleID: A Column of sampleIDs use for Training
# chr: Chromosome Number for Corresponding Genotype Data
# genofile_type: vcf or dosages
# genofile: Genotype File path, should be tabix (contains .gz and .tbi)
# Format: Format using for Genotype Data (GT or DS).
# maf: Threshold for Minor Allele Frequency (range from 0-1). Default is 0.01.
# hwe: Threshold of p-value for Hardy Weinberg Equilibrium exact test. Default is 0.001.
# window: Window size around gene boundary. Default is 10**6 BP.
# thread: Number of Thread for multiprocessing
# out: Output Dir

### Elastic Net
# cv: cv-fold cross-validation in model selection, default 5-fold
# alpha: L1 & L2 ratio for elastic net regression, default 0.5
#        If alpha=0, lasso regression
#        If alpha=1, ridge regression

### DPR
# dpr: model using for DPR
# ES: Effect Size (fixed or additive). Default is fixed.

#############################################################################################################
VARS=`getopt -o "" -a -l \
model:,Gene_Exp:,sampleID:,chr:,genofile_type:,genofile:,Format:,maf:,hwe:,window:,cv:,alpha:,dpr:,ES:,thread:,out: \
-- "$@"`

if [ $? != 0 ]
then
    echo "Terminating....." >&2
    exit 1
fi
 
eval set -- "$VARS"

while true
do
    case "$1" in
        --model|-model) model=$2; shift 2;;
        --Gene_Exp|-Gene_Exp) Gene_Exp=$2; shift 2;;
        --sampleID|-sampleID) sampleID=$2; shift 2;;
        --chr|-chr) chr_num=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --Format|-Format) Format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --hwe|-hwe) hwe=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --cv|-cv) cv=$2; shift 2;;
        --alpha|-alpha) alpha=$2; shift 2;;
        --dpr|-dpr) dpr_num=$2; shift 2;;
        --ES|-ES) ES=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### Setting Default Values
maf=${maf:-0.01}
hwe=${hwe:-0.001}
window=${window:-$((10**6))}
thread=${thread:-1}

# model==elastic_net
cv=${cv:-5}
alpha=${alpha:-0.5}

# model==DPR
dpr_num=${dpr_num:-1}
ES=${ES:-"fixed"}

##########################################################################################################
### Checking --model Input
if [[ "$model"x != "elastic_net"x ]] && [[ "$model"x != "DPR"x ]]; then
    echo "Model Not Found."
    echo "Please input elastic_net or DPR for --model option."
else 
    ######################################################################################################
    ### 1. 
    ### Create Dir & Store Result
    mkdir -p ${out_prefix}/${model}_CHR${chr_num}

    ### Store python log files
    mkdir -p ${out_prefix}/${model}_CHR${chr_num}/log_file

    ### 2. 
    ### Extract Genotype File Header
    zcat ${genofile} | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt

    ### 3.
    ### Extract Gene and Expression Level Data by Chromosome
    grep 'CHROM' ${Gene_Exp} > ${out_prefix}/CHR${chr_num}_Gene_Exp_combination.txt
    awk '{if ($1=="'"$chr_num"'") print}' ${Gene_Exp} >> ${out_prefix}/CHR${chr_num}_Gene_Exp_combination.txt

    ### 4.
    ### Training
    ##################################################################################################
    ### Elastic Net Regression
    ##################################################################################################
    if [[ "$model"x == "elastic_net"x ]];then
        echo "Using Elastic Net Model for Training."

        ### 3.1 Elastic Net Regression
        python ./Model_Train_Pred/Elastic_Net_Train.py \
        --Gene_Exp_path ${out_prefix}/CHR${chr_num}_Gene_Exp_combination.txt \
        --sampleID ${sampleID} \
        --chr_num ${chr_num} \
        --genofile ${genofile} \
        --genofile_header ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt \
        --geno ${genofile_type} \
        --Format ${Format} \
        --maf ${maf} \
        --hwe ${hwe} \
        --window ${window} \
        --cv ${cv} \
        --alpha ${alpha} \
        --thread ${thread} \
        --out_prefix ${out_prefix}/${model}_CHR${chr_num} \
        > ${out_prefix}/${model}_CHR${chr_num}/log_file/${model}_Train_log.txt

    ##################################################################################################
    ### DPR
    ##################################################################################################

    elif [[ "$model"x == "DPR"x ]]; then
        echo "Using DPR Model for Training."

        ### 3.3
        ### Store DPR Input Files
        mkdir -p ${out_prefix}/${model}_CHR${chr_num}/DPR_input

        ### Store bimbam file
        mkdir -p ${out_prefix}/${model}_CHR${chr_num}/DPR_input/bimbam

        ### Store phenotype file
        mkdir -p ${out_prefix}/${model}_CHR${chr_num}/DPR_input/pheno

        ### Store SNP annotation file
        mkdir -p ${out_prefix}/${model}_CHR${chr_num}/DPR_input/SNP_annot

        ### Store cross validation result
        mkdir -p ${out_prefix}/${model}_CHR${chr_num}/DPR_input/CV 
        mkdir -p ${out_prefix}/${model}_CHR${chr_num}/DPR_input/CV/bimbam
        mkdir -p ${out_prefix}/${model}_CHR${chr_num}/DPR_input/CV/pheno
        mkdir -p ${out_prefix}/${model}_CHR${chr_num}/DPR_input/CV/SNP_annot

        ### 3.4 DPR
        python ./Model_Train_Pred/DPR_Train.py \
        --Gene_Exp_path ${out_prefix}/CHR${chr_num}_Gene_Exp_combination.txt \
        --sampleID ${sampleID} \
        --chr_num ${chr_num} \
        --genofile ${genofile} \
        --genofile_header ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt \
        --geno ${genofile_type} \
        --Format ${Format} \
        --hwe ${hwe} \
        --maf ${maf} \
        --window ${window} \
        --dpr ${dpr_num} \
        --ES ${ES} \
        --thread ${thread} \
        --out_prefix ${out_prefix}/${model}_CHR${chr_num} \
        > ${out_prefix}/DPR_CHR${chr_num}/log_file/${model}_Train_log.txt

        ### 3.5 Remove Files
        rm -r ${out_prefix}/${model}_CHR${chr_num}/DPR_input/CV
##################################################################################################
    fi
fi

##################################################################################################
### 5.
### Remove Files
rm ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt
rm ${out_prefix}/CHR${chr_num}_Gene_Exp_combination.txt











