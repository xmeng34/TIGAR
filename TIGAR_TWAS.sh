#!/bin/bash

#############################################################################################################
# software requirement
# python 3
# tabix
###

###
# asso: Method for association study:
#       1) asso==1, run TWAS with PED file
#       2) asso==2, run TWAS with provided GWAS Z-score

# Gene_Exp: Predicted gene expression data
# thread : Number of Thread. Default is 1.
# out : Path for TIGAR to save output files.

# For asso==1
# PED : PED File.
# Asso_Info : Instruction for association study
#             1) P : columns' names for corresponding phenotype in PED file
#             2) C : columns' names for covariates in PED file
# method : Link Function. Default is OLS.
#          1) OLS : Ordinary Least Square Regression
#          2) Logit : Logistic Regression

# For asso==2
# Zscore : Zscore file from previous GWAS study (tabixed)
# Weight : File contains snps effect size (Same Format as Training Output File).
# Covar : Reference Covariance Matrix (Scripts is provided, see in covar_calculation.py, TIGAR_Covar.sh, tabixed).
# chr : Chromosome Number
# window : Window Size around Gene Boundary. Default is 10^6BP.

#############################################################################################################
VARS=`getopt -o "" -a -l \
asso:,Gene_Exp:,PED:,Asso_Info:,method:,Zscore:,Weight:,Covar:,chr:,window:,thread:,out: \
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
        --asso|-asso) asso=$2; shift 2;;
        --Gene_Exp|-Gene_Exp) Gene_Exp=$2; shift 2;;
        --PED|-PED) PED=$2; shift 2;;
        --Asso_Info|-Asso_Info) Asso_Info=$2; shift 2;;
        --method|-method) method=$2; shift 2;;
        --Zscore|-Zscore) Zscore=$2; shift 2;;
        --Weight|-Weight) Weight=$2; shift 2;;
        --Covar|-Covar) Covar=$2; shift 2;;
        --chr|-chr) chr_num=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### s=Setting Default Values
thread=${thread:-1}

### asso==1
method=${method:-'OLS'}
PED=${PED:-""}
Asso_Info=${Asso_Info:-""}

### asso==2
window=${window:-$((10**6))}
Zscore=${Zscore:-""}
Weight=${Weight:-""}
Covar=${Covar:-""}
chr_num=${chr_num:-""}

#####################################################################################
# TWAS
#####################################################################################
# Checking Association Type Input (--asso)
if [[ "$asso" -ne 1 ]] && [[ "$asso" -ne 2 ]]; then
    echo "Type of Association Study Not Found."
    echo "Please input 1 or 2 for --asso"
else 
   #####################################################################################
   # Type One Association Study
   #####################################################################################
    if [[ "$asso" -eq 1 ]]; then
        echo "Type One Association Study (asso is 1)"

        mkdir -p ${out_prefix}/TIGAR_TWAS_Type_One
        mkdir -p ${out_prefix}/TIGAR_TWAS_Type_One/log_file
        
        python ./TWAS/TWAS.py ONE \
        --Gene_Exp_Path ${Gene_Exp} \
        --PED ${PED} \
        --Asso_Info ${Asso_Info} \
        --method ${method} \
        --thread ${thread} \
        --out_prefix ${out_prefix}/TIGAR_TWAS_Type_One \
        > ${out_prefix}/TIGAR_TWAS_Type_One/log_file/${method}_PHENO_NUM_`cut -f1 ${Asso_Info} | grep 'P' | wc -l`_Type_One.txt

    ####################################################################################
    # Type Two Association Study
    ####################################################################################
    elif [[ "$asso" -eq 2 ]]; then
        echo "Type Two Association Study (asso is 2)"

        mkdir -p ${out_prefix}/TIGAR_TWAS_Type_Two
        mkdir -p ${out_prefix}/TIGAR_TWAS_Type_Two/log_file

        ### Tabix Training Weight File
        sed -n '2,$p' ${Weight} | sort -n -k2 | bgzip -c > ${out_prefix}/CHR${chr_num}_weight.txt.gz
        tabix -p vcf ${out_prefix}/CHR${chr_num}_weight.txt.gz
        ### Extract Header of Weight File
        grep 'CHROM' ${Weight} > ${out_prefix}/CHR${chr_num}_weight_header.txt

        ### Extract Gene Annotation File
        echo -e 'CHROM\tGeneStart\tGeneEnd\tTargetID\tGeneName' > ${out_prefix}/CHR${chr_num}_Gene_Annotation.txt
        awk '{if ($1=="'"$chr_num"'") print}' ${Gene_Exp} | cut -f1-5 \
        >> ${out_prefix}/CHR${chr_num}_Gene_Annotation.txt

        ### Extract Header of Z-score File
        zcat ${Zscore} | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_Zscore_header.txt

        python ./TWAS/TWAS.py TWO \
        --Gene ${out_prefix}/CHR${chr_num}_Gene_Annotation.txt \
        --Zscore ${Zscore} \
        --Zscore_header ${out_prefix}/CHR${chr_num}_Zscore_header.txt \
        --Weight ${out_prefix}/CHR${chr_num}_weight.txt.gz \
        --Weight_header ${out_prefix}/CHR${chr_num}_weight_header.txt \
        --Covar ${Covar} \
        --chr_num ${chr_num} \
        --window ${window} \
        --thread ${thread} \
        --out_prefix ${out_prefix}/TIGAR_TWAS_Type_Two \
        > ${out_prefix}/TIGAR_TWAS_Type_Two/log_file/CHR${chr_num}_Type_Two_log.txt

        ### Remove Files
        rm ${out_prefix}/CHR${chr_num}_Gene_Annotation.txt
        rm ${out_prefix}/CHR${chr_num}_weight_header.txt
        rm ${out_prefix}/CHR${chr_num}_weight.txt.gz
        rm ${out_prefix}/CHR${chr_num}_weight.txt.gz.tbi
        rm ${out_prefix}/CHR${chr_num}_Zscore_header.txt
######################################################################################
    fi
fi





































