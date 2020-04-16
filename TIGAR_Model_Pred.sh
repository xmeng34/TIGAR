#!/bin/bash

#####################################################################################################
# software requirement
# python 3
# tabix
###

### Variable Needed for Prediction
# model: elastic_net or DPR
# chr: Chromosome Number for Corresponding Genotype Data
# train_weight_path: Training Weight File (Same Format as Training Output)
# train_info_path: Training Information File (Same Format as Training Output)
# genofile_type: vcf or dosages
# genofile: Genotype File path. Should be tabixed. (contains .gz and .tbi)
# sampleID: A Column of sampleIDs use for Prediction.
# Format: Format using for Genotype Data (GT or DS).
# window: Window size around gene boundary. Default is 10**6 BP.
# maf_diff: Threshold of difference between maf calculated in Training and Prediction. 
#           If difference is larger than this value, then drop the snp. Default is 0.2. 
# thread: Number of Thread for multiprocessing
# out: Output Dir

#####################################################################################################
VARS=`getopt -o "" -a -l \
model:,chr:,train_weight_path:,train_info_path:,genofile_type:,genofile:,sampleID:,Format:,window:,maf_diff:,thread:,out: \
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
        --chr|-chr) chr_num=$2; shift 2;;
        --train_weight_path) train_weight_path=$2; shift 2;;
        --train_info_path|-train_info_path) train_info_path=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --sampleID|-sampleID) sampleID=$2; shift 2;;
        --Format|-Format) Format=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --maf_diff|-maf_diff) maf_diff=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### Setting Default Values
window=${window:-$((10**6))}
maf_diff=${maf_diff:-0.2}
thread=${thread:-1}

###############################################################################
# Checking --model input
if [[ "$model"x != "elastic_net"x ]] && [[ "$model"x != "DPR"x ]]; then
    echo "Model Not Found."
    echo "Please input elastic_net or DPR for --model."
else 
    ### Extract Genotype File Header
    zcat ${genofile} | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt
    
    ### Store Results
    mkdir -p ${out_prefix}/${model}_CHR${chr_num}
    mkdir -p ${out_prefix}/${model}_CHR${chr_num}/log_file
    
    ### Prediction
    python ./Model_Train_Pred/Prediction.py \
    --model ${model} \
    --chr_num ${chr_num} \
    --train_weight_path ${train_weight_path} \
    --train_info_path ${train_info_path} \
    --genofile ${genofile} \
    --genofile_header ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt \
    --sampleID ${sampleID} \
    --geno ${genofile_type} \
    --Format ${Format} \
    --window ${window} \
    --maf_diff ${maf_diff} \
    --thread ${thread} \
    --out_prefix ${out_prefix}/${model}_CHR${chr_num} \
    > ${out_prefix}/${model}_CHR${chr_num}/log_file/${model}_Prediction_log.txt
    
    ### Remove Files
    rm ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt

fi







