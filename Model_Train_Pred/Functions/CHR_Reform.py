#!/usr/bin/env python

import numpy as np 
from numpy import *

from dfply import *

import sys
sys.path.append("./Model_Train_Pred/Functions")
from HWE import HWE

#######################################################################################################
### Input:
### 1. part:
###    1) TRAIN: Used in Model Training Part
###    2) PRED: Used in Model Prediction Part
### 2. data: Genotype Data Extracted by Tabix
### 3. geno: vcf or dosages
### 4. sampleID: sampleID used in Training or Prediction. Based on Genotype Data Provided.
### 5. Format: GT or DS
### 6. p_hwe: Threshold of p-value for Hardy Weinberg Equilibrium Exact Test
### 7. maf: Threshold of Minor Allele Frequency (Range from 0-1)

### Output:
### 1. The First Six Columns of Output Data Frame should be:
###    1) CHROM
###    2) POS
###    3) ID (i.e. rsID)
###    4) REF
###    5) ALT
###    6) snpID (CHROM:POS:REF:ALT)
###    7) p_HWE: p-value for Hardy Weinberg Equilibrium Exact Test
###    8) MAF: Minor Allele Frequency (Range from 0-1)
###    9) Samples' Genotype Data Split by Format (GT or DS)

### 2. Output Data Frame will be Selected by p_HWE and MAF:
###    1) TRAIN: Default threshold is p_HWE > 10**(-3) and MAF > 0.01
###    2) PRED: MAF > 0

####################################################################################################
class CHR_Reform(object):
    def __init__(self, part, data, geno, sampleID, Format, p_hwe=-1, maf=-1):
        self.part = part
        self.data = data
        self.geno = geno
        self.sampleID = sampleID
        self.Format = Format
        self.p_hwe = p_hwe
        self.maf = maf

    ################################################################################
    ### Input Genotype Data by SampleID
    ### For GT Format:
    ###  code '0|0' or '0/0' as 0
    ###  code ('0|1' or '1|0')  or ('0/1' or '1/0') as 1
    ###  code '1|1' or '1/1' as 2
    ###  code '.|.' or './.' as nan (missing)

    ### For DS Format:
    ### code '.' as nan (missing)
    ##################################################################################
    def geno_reform(self, data, Format):
        if Format == 'GT':
            data[(data=='0|0')|(data=='0/0')]=0
            data[(data=='1|0')|(data=='1/0')|(data=='0|1')|(data=='0/1')]=1
            data[(data=='1|1')|(data=='1/1')]=2
            data[(data=='.|.')|(data=='./.')]=nan
        elif Format == 'DS':
            data[(data=='.')]=nan
        return data

    def main(self):
        if self.geno == 'vcf':
            sampleID = np.intersect1d(self.data.columns[9:], self.sampleID)
        elif self.geno == 'dosages':
            sampleID = np.intersect1d(self.data.columns[5:], self.sampleID)

        ### Generate snpIDs
        self.data['snpID']=(self.data['CHROM'].astype('str')+":"+self.data['POS'].astype('str')
                            +":"+self.data.REF+":"+self.data.ALT)

        CHR = self.data >> select(self.data[['CHROM','POS','ID','REF','ALT','snpID']],
                                  self.data[sampleID])

        CHR=CHR.drop_duplicates(['snpID'],keep='first')

        if self.geno == 'vcf':
            indicate=self.data.FORMAT[0].split(":").index(self.Format)
            CHR[sampleID]=CHR[sampleID].applymap(lambda x:x.split(":")[indicate])

        CHR[sampleID]=CHR[sampleID].apply(lambda x:self.geno_reform(x,self.Format),axis=0)

        ### Calculate p_HWE (Optional) and MAF by SNPs
        temp=pd.DataFrame((CHR >> select(CHR[sampleID])),dtype=np.float)

        ### Calculating p_HWE
        if self.part == 'TRAIN':
            from HWE import HWE
            CHR['p_HWE']=temp.apply(lambda x:HWE(x.dropna()).main(),axis=1)

        ### Calculate MAF(Range from 0-1)
        CHR['MAF']=temp.apply(lambda x:sum(x)/(2*len(x.dropna())),axis=1)
        ### Dealing with NaN
        CHR[np.hstack(([sampleID,'MAF']))] = CHR[np.hstack(([sampleID,
                                                             'MAF']))].apply(lambda x:x.fillna(2*x.MAF),axis=1)

        if self.part == 'TRAIN':
            return (CHR >> mask((CHR.p_HWE>=self.p_hwe)&(CHR.MAF>=self.maf))).reset_index(drop=True)
        elif self.part == 'PRED':
            CHR = CHR >> mask(CHR.MAF>=0)
            return (CHR.rename(columns={'MAF':'MAF_test'})).reset_index(drop=True)












