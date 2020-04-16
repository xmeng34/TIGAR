#!/usr/bin/env python

####################################################################
### Import Package Needed
import argparse
import warnings
from time import *

import pandas as pd
import numpy as np
from numpy import *
from dfply import *

import sys
sys.path.append("./TWAS/Functions")

################################################################################
# Hide Warnings from Python
warnings.filterwarnings("ignore")
################################################################################
### Time Calculation
start_time=time()

################################################################################
# Type One Association Study
#################################################################################
def Type_One(args):
    print("Type One Association Study.")

    ### Checking Input Commands
    print("Gene Annotation and Expression Level File:"+args.Gene_Exp_Path)
    print("PED File:"+args.PED)
    print("Association Information File:"+args.Asso_Info)
    print("Link Functions:"+args.method)
    print("Number of Thread:"+str(args.thread))
    print("Output Dir:"+args.out_prefix)

    # Read in PED file
    PED = pd.read_csv(args.PED,sep='\t').rename(columns={'#FAM_ID':'FAM_ID'})

    # Gene Annotation and Expression Level File
    Gene_Exp = pd.read_csv(args.Gene_Exp_Path,sep='\t')

    # Read in Association information
    # P : Phenotype
    # C : Covariate
    Asso_Info=pd.read_csv(args.Asso_Info,sep='\t',header=None)
    Asso_Info.columns=['Ind','Var']

    ### Phenotype
    pheno = Asso_Info >> mask(Asso_Info.Ind=='P') >> select(Asso_Info.Var) 
    pheno = np.array(pheno).ravel()
    if len(pheno)==0:
        raise SystemExit("No Phenotype Provided.")
    else:
        print("Phenotype Used:",pheno)

    ### Covariances
    covar = Asso_Info >> mask(Asso_Info.Ind=='C') >> select(Asso_Info.Var)
    covar = np.array(covar).ravel()
    if len(covar)==0:
        raise SystemExit("No Covariates Provided.")
    else:
        print("Covariance Used:",covar)
    
    ### TargetID
    TargetID = np.array(Gene_Exp.TargetID)
    if len(TargetID)==0:
        raise SystemExit("No Gene in Gene-Expression File.")
    
    ### sampleIDs
    sampleID = np.intersect1d(np.array(PED.IND_ID),np.array(Gene_Exp.columns[5:]))
    if len(sampleID)==0:
        raise SystemExit("No Overlapping IND_ID between Gene-Expression and PED File.")

    # Organizing PED and Gene-Expression File
    PED = PED >> mask(PED.IND_ID.isin(sampleID)) >> select(PED.IND_ID,PED[pheno],PED[covar])
    # Code 'X' in PED File with nan (missing)
    PED[PED=='X']=np.nan
    Gene_Exp=Gene_Exp >> select(Gene_Exp[Gene_Exp.columns[0:5]],Gene_Exp[sampleID])

    if len(pheno) == 1:
        from Single_Phenotype import Single_Phenotype
        print("Single Phenotype.")

        Single_Phenotype(args.method, pheno, covar, PED, Gene_Exp, 
                         TargetID, args.out_prefix, args.thread).main()
    elif len(pheno) > 1:
        from Multi_Phenotype import Multi_Phenotype
        print("Multiple Phenotype.")

        Multi_Phenotype(args.method, pheno, covar, PED, Gene_Exp,
                        TargetID, args.out_prefix, args.thread).main()

################################################################################
# Type Two Association Study
#################################################################################
def Type_Two(args):
    print("Type Two Association Study.")

    ### Checking Input Commands
    print("Gene Annotation File:"+args.Gene)
    print("Z-score from Previous GWAS Study:"+args.Zscore)
    print(args.Zscore_header)
    print("Snps Weight:"+args.Weight)
    print(args.Weight_header)
    print("Reference Covariance File:"+args.Covar)
    print("Chromosome Number:"+str(args.chr_num))
    print("window="+str(args.window))
    print("Number of Thread:"+str(args.thread))
    print("Output Dir:"+args.out_prefix)
    
    ### Read in Gene Annotation 
    Gene = pd.read_csv(args.Gene,sep='\t')

    ### TargetID
    TargetID = np.array(Gene.TargetID)

    ### Initializing Header of Output File
    pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd', 'TargetID',
                          'GeneName','ZSCORE','PVALUE']).to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_association_study.txt',
                                                                sep='\t',index=None,header=True,mode='a')
    
    ### Read in Header for Weight & Zscore Files
    Weight_header = pd.read_csv(args.Weight_header,sep='\t')
    Zscore_header = pd.read_csv(args.Zscore_header,sep='\t')
    
    from Summary_Stat import Summary_Stat
    Summary_Stat(Gene, args.Zscore, Zscore_header, args.Weight, Weight_header, args.Covar,
                 args.chr_num, TargetID, args.window, args.out_prefix, args.thread).main()

################################################################################
# Inputs
#################################################################################
parser = argparse.ArgumentParser(prog='Association Study')

subparsers = parser.add_subparsers(help='sub-command help')
#################################################################################
### ONE
### variables need for Type One Association Study
parser_one = subparsers.add_parser('ONE', help='Type One Association Study')

### Gene Annotation and Expression Level File
parser_one.add_argument('--Gene_Exp_Path', type=str, default=None)

### PED File
parser_one.add_argument('--PED', type=str, default=None)

### Association Information File (P & C columns)
parser_one.add_argument('--Asso_Info', type=str, default=None)

### Link Function (OLS or Logit)
parser_one.add_argument('--method', type=str, choices=['OLS', 'Logit'], default=None, 
                        help='Link Function.')

### Number of Thread
parser_one.add_argument('--thread', type=int, default=None)

### Output Dir
parser_one.add_argument('--out_prefix', type=str, default=None)

parser_one.set_defaults(func=Type_One)

#################################################################################
### TWO
### variables need for Type Two Association Study
parser_two = subparsers.add_parser('TWO', help='Type Two Association Study')

### Gene Annotation File
parser_two.add_argument('--Gene', type=str, default=None)

###  Zscore File from Previous GWAS Study (tabixed)
parser_two.add_argument('--Zscore', type=str, default=None)
parser_two.add_argument('--Zscore_header', type=str, default=None)

### File Contains snps Effect-Size (Same Format as Training Output File, tabixed)
parser_two.add_argument('--Weight', type=str, default=None)
parser_two.add_argument('--Weight_header', type=str, default=None)

### Reference Covariance File (tabixed)
parser_two.add_argument('--Covar', type=str, default=None)

### Chromosome Number
parser_two.add_argument('--chr_num', type=int, default=None)

### Window Size around Gene Boundary
parser_two.add_argument('--window', type=int, default=None)

### Number of Thread
parser_two.add_argument('--thread', type=int, default=None)

### Output Dir
parser_two.add_argument('--out_prefix', type=str,default=None)

parser_two.set_defaults(func=Type_Two)

#################################################################################
args = parser.parse_args()

### Execute Function
args.func(args)

#################################################################################
### Time Calculation
end_time=time()
print("Running Time:"+str(round((end_time-start_time)/60,2))+" minutes.")





