## TIGAR
"TIGAR" standing for Transcriptome-Integrated Genetic Association Resource, which is developed using Python and BASH scripts. TIGAR can fit both Elastic-Net and nonparametric Bayesian model (Dirichlet Process Regression, i.e. DPR) for gene expression imputation, impute genetically regulated gene expression (GReX) from genotype data, and conduct transcriptome-wide association studies (TWAS) using both individual-level and summary-level GWAS data for univariate and multivariate phenotypes.

### Software
1. DPR
	- DPR module is saved under folder `./Model_Train_Pred/Functions`

	- Please add the executable file `./Model_Train_Pred/DPR` to your linux `${PATH}` directory. Assuming `~/bin/` is a directory added to your `${PATH}` environmental variable, you can accommodate the following example command

		```
		cp ./Model_Train_Pred/Functions/DPR ~/bin/
		```
	OR

	- Please run the following command before running DPR in TIGAR.
		```
		chmod a+x ./Model_Train_Pred/Functions/DPR
		```

2. BGZIP: http://www.htslib.org/doc/bgzip.html 
3. TABIX: http://www.htslib.org/doc/tabix.html 
4. python 3.5 
   - dfply
   - io
   - subprocess
   - multiprocess


### Input File Format
Example data provided here are generated artificially. All input files are tab delimited text files.


#### 1. Gene Annotation and Expression File (`./example_data/Gene_Exp.txt`)
- First 5 columns specify chromosome number, gene start position, gene end position, target gene ID, gene name (optional, could be the same as gene ID).
- Sample gene expression data start from the 6th column. 

| CHROM | GeneStart | GeneEnd |   TargetID      | GeneName | sample1 | sample...|
|:-----:|:---------:|:-------:|:---------------:|:--------:|:-------:|:--------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |   0.2   |     ...  |


#### 2. Genotype File
- Sorted by chromosome and base pair position, zipped by `bgzip`, and tabixed.
- Example tabix command, `tabix -f -p vcf *.vcf.gz`.

- vcf file (`./example_data/Genotype/example.vcf.gz`)
	- Genotype data start from the 10th column.
	- More information about VCF file format: http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/

	| CHROM | POS |  ID | REF | ALT | QUAL | FILTER | INFO | FORMAT |  sample1 | sample...|
	|:-----:|:---:|:---:|:---:|:---:|:----:|:------:|:----:|:------:|:--------:|:--------:|
	|   1   | 100 | rs1 |  C  |  T  |   .  |  PASS  |   .  |  GT:DS | 0/0:0.01 |    ...   |

- dosages file
	- The first 5 columns are of the same format as VCF file.
	- Dosage genotype data start from the 6th column.

	| CHROM | POS |  ID | REF | ALT | sample1 | sample...|
	|:-----:|:---:|:---:|:---:|:---:|:-------:|:--------:|
	|   1   | 100 | rs1 |  C  |  T  |   0.01  |    ...   |

#### 3. PED File (`./example_data/example_PED.ped`)
- More information bout PED file format: http://zzz.bwh.harvard.edu/plink/data.shtml#ped

| FAM_ID | IND_ID | FAT_ID | MOT_ID | SEX | PHENO | COV1 | COV...|
|:------:|:------:|:------:|:------:|:---:|:-----:|:---:|:---:|
|   11A  |   11A  |    X   |    X   |  1  |  0.2  | 0.3 |...|

#### 4. Asso_Info File (`./example_data/Asso_Info/Asso_Info_*.txt`)
- Two columns with the first column specifying the Phenotype (P) and Covariate variables (C) from the PED file, and the second column specifying the corresponding variable names in the PED file. The variables specified in the Asso_Info file will be used in TWAS.

|P|PHENO|
|:-----:|:---:|
|C|COV1|
|C|COV2|
|C|SEX|

#### 5.Zscore File (`./example_data/example_Zscore/*_GWAS_Zscore.txt.gz`)
- Sorted by chromosome and base pair position, zipped by `bgzip`, and tabixed
- Example tabix command, `tabix -f -p vcf *_Zscore.txt.gz`. The first 4 columns are of the same format as VCF file.

| CHROM | POS | REF | ALT | Zscore |
|:-----:|:---:|:---:|:---:|:------:|
|   1   | 100 |  C  |  T  |  0.01  |

#### 6. Weight File used for TWAS with GWAS summary statistics
- First 5 columns have to be of the following format, specifying chromosome number, base pair position, reference allele, alternative allele, and target gene ID. 

- The column `ES` (Effect Size) denotes the weights for this given SNP/Target(Gene)

| CHROM | POS | REF | ALT |     TargetID    |  ES  |
|:-----:|:---:|:---:|:---:|:---------------:|:----:|
|   1   | 100 |  C  |  T  |     ENSG0000    |  0.2 |

#### 7. Genome Block Annotation File (`./example_data/block_annotation_EUR.txt`)
- The block annotation file is a tab delimited text file with head row of `CHROM Start End`, denoting the chromosome number, starting position, ending position. Reference genotype files shall be of one per chromosome, or one for the whole genome-wide variants. Example block annotation file for European samples is provided `./TIGAR/example_data/block_annotation_EUR.txt`. 

| CHROM |   Start   |   End   |
|:-----:|:---------:|:-------:|
|   1   |    100    | 20000   |  

- Block annotation files of other ethnicities can be adopted from the genome segmentation generated by `LDetect`, https://bitbucket.org/nygcresearch/ldetect-data/src/master/.

### Example Usage 
1. More details are available in the TIGAR_Manual.pdf

2. cis-eQTL effect-size calculation

- Train Elastic-Net Imputation Model
```
Gene_Exp_path=./example_data/Gene_Exp.txt
sampleID=./example_data/sampleID.txt
genofile=./example_data/Genotype/example.vcf.gz
out_prefix=./Result

cd TIGAR

./TIGAR_Model_Train.sh --model elastic_net \
--Gene_Exp ${Gene_Exp_path} --sampleID ${sampleID} \
--chr 1 --genofile_type vcf \
--genofile ${genofile} --Format GT \
--out ${out_prefix}
```

3. GReX Prediction
```
genofile=./example_data/Genotype/example.vcf.gz
sampleID=./example_data/sampleID.txt
train_weight_path=./Result/elastic_net_CHR1/CHR1_elastic_net_training_weight.txt
train_info_path=./Result/elastic_net_CHR1/CHR1_elastic_net_training_info.txt
out_prefix=./Result

cd TIGAR

./TIGAR_Model_Pred.sh --model elastic_net \
--chr 1 \
--train_weight_path ${train_weight_path} \
--train_info_ptah ${train_info_path} \
--genofile_type vcf \
--genofile ${genofile} --Format GT \
--out ${out_prefix}
```

4. TWAS
- Using individual-level GWAS data. Take the output `*_GReX_prediction.txt` from gene expression prediction as the input for `--Gene_EXP` here.

```
Gene_Exp_path=./Result/DPR_CHR1/CHR1_DPR_GReX_prediction.txt
PED=./example_data/example_PED.ped
Asso_Info=./example_data/Asso_Info/Asso_Info_SinglePheno_OLS.txt
out_prefix=./Result/DPR_CHR1

cd TIGAR

./TIGAR_TWAS.sh --asso 1 \
--Gene_EXP ${Gene_Exp_path} \
--PED ${PED} --Asso_Info ${Asso_Info} \
--out ${out_prefix}
```

- Using summary-level GWAS data. Take the output `*_training_weight.txt` from imputation model training as the input Weight file here. Take original gene annotation and expression level (`./example_data/Gene_Exp.txt`) as input for `--Gene_Exp`

```
Gene_Exp_path=./example_data/Gene_Exp.txt
Zscore=./example_data/example_Zscore/CHR1_GWAS_Zscore.txt.gz
Weight=./DPR_CHR1/CHR1_DPR_training_param.txt
Covar=./Result/reference_cov/CHR1_reference_cov.txt.gz
out_prefix=./Result/DPR_CHR1

cd TIGAR

./TIGAR_TWAS.sh --asso 2 \
--Gene_Exp ${Gene_Exp_path} 
--Zscore ${Zscore} --Weight ${Weight} --Covar ${Covar} --chr 1 \
--out ${out_prefix}
```

- Generate Reference Covariance File
```
block=./example_data/block_annotation_EUR.txt
genofile=./example_data/example.vcf.gz
out_prefix=./Result

cd TIGAR

./TWAS/Covar/TIGAR_Covar.sh --block ${block} \
--genofile_type vcf --genofile ${genofile} \
--chr 1 \
--Format GT \
--out ${out_prefix}
```

### Reference
- Elastic Net: https://github.com/hakyimlab/PrediXcan  
- DPR: https://github.com/biostatpzeng/DPR






