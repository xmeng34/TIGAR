#!/bin/bash

#####################################################################################################
# software requirement
# python 3
# tabix
###

###
# block : Block Annotation File
# genofile_type : vcf or dosages
# genofile: Genofile path, should be tabix(contains .gz and .tbi)
# chr : Chromosome Number
# Format: Format Using for Genotype Data (GT or DS).
# maf : Threshold for Minor Allele Frequency (range from 0-1),default 0.05
# thread : Number of thread
# out : output dir
#####################################################################################################

VARS=`getopt -o "" -a -l \
block:,genofile_type:,genofile:,chr:,Format:,maf:,thread:,out: \
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
        --block|-block) block=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --chr|-chr) chr_num=$2; shift 2;;
        --Format|-Format) Format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### setting default value
thread=${thread:-1}
maf=${maf:-0.05}

###############################################################################################
### 1. Store Results
mkdir -p ${out_prefix}/reference_cov

### 2. Extract genofile header
zcat ${genofile} | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt

### 3. Extract Block Annotation by Chromosome
echo -e 'CHROM\tStart\tEnd' > ${out_prefix}/CHR${chr_num}_block_annotation.txt
awk '{if ($1=="'"$chr_num"'") print}' ${block} >> ${out_prefix}/CHR${chr_num}_block_annotation.txt

### 4. Run covariance calculation
python ./TWAS/Covar/covar_calculation.py \
--block ${out_prefix}/CHR${chr_num}_block_annotation.txt \
--genofile ${genofile} \
--genofile_header ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt \
--geno ${genofile_type} \
--chr_num ${chr_num} \
--Format ${Format} \
--maf ${maf} \
--thread ${thread} \
--out ${out_prefix}/reference_cov \
> ${out_prefix}/reference_cov/CHR${chr_num}_Covar_log.txt

### 2. tabix output file
sort -n -k2 ${out_prefix}/reference_cov/CHR${chr_num}_reference_cov.txt | bgzip -c \
> ${out_prefix}/reference_cov/CHR${chr_num}_reference_cov.txt.gz

tabix -p vcf ${out_prefix}/reference_cov/CHR${chr_num}_reference_cov.txt.gz

### 3. Remove Files
rm ${out_prefix}/CHR${chr_num}_${genofile_type}_header.txt
rm ${out_prefix}/reference_cov/CHR${chr_num}_reference_cov.txt
rm ${out_prefix}/CHR${chr_num}_block_annotation.txt














