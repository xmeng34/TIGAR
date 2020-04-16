#!/usr/bin/env python

################################################################################
# Import Packages Need
import subprocess
from subprocess import *
import pandas as pd
import numpy as np
from numpy import *
from dfply import *
import io
from io import StringIO
from io import *
from scipy.stats import chi2

import multiprocessing

################################################################################
### Input:
### 1. Gene: Gene Annotation File
### 2. Zscore: Path of Z-score File 
### 3. Zscore_header: Z-score Header
### 4. Weight: Path of Weight File
### 5. Weight_header: Weight Header
### 6. Covar: Reference Covariance File (tabixed)
### 7. chr_num: Chromosome Number
### 8. TargetID
### 9. window: Window Size around Gene Boundary
### 10. out_prefix: Output Dir
### 11. thread: Number of Thread

### Output:
### 1. CHROM, GeneStart, GeneEnd ,GeneName, TargetID: Gene Annotation columns
### 2. ZSCORE: Value of Burden Z-score
### 3. PVALUE: P-value of Chi-Square Test for Burden Z-score

################################################################################
class Summary_Stat(object):
    def __init__(self, Gene, Zscore, Zscore_header, Weight, Weight_header, Covar, 
                 chr_num, TargetID, window, out_prefix, thread):
        self.Gene = Gene
        self.Zscore = Zscore
        self.Zscore_header = Zscore_header
        self.Weight = Weight
        self.Weight_header = Weight_header
        self.Covar = Covar
        self.chr_num = chr_num
        self.TargetID = TargetID
        self.window = window
        self.out_prefix = out_prefix
        self.thread = thread
        
    def thread_process(self, num):
        Gene_temp = self.Gene >> mask(self.Gene.TargetID == self.TargetID[num])
        
        start = max(int(Gene_temp.GeneStart) - self.window, 0)
        end = max(int(Gene_temp.GeneEnd) + self.window, 0)
        
        Zscore_process = subprocess.Popen(["tabix"+" "+self.Zscore+" "+str(self.chr_num)+":"+str(start)+"-"+str(end)],
                                          shell=True,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
        Zscore_out = Zscore_process.communicate()[0]

        Weight_process = subprocess.Popen(["tabix"+" "+self.Weight+" "+str(self.chr_num)+":"+str(start)+"-"+str(end)],
                                          shell=True,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
        Weight_out = Weight_process.communicate()[0]
        
        if len(Zscore_out)*len(Weight_out) == 0:
            print("No GWAS Test Z-score or Weight for this Gene:"+self.TargetID[num])
        
        else:
            print("Running Association Study for Gene:"+self.TargetID[num])
            Zscore = pd.read_csv(StringIO(Zscore_out.decode('utf-8')),sep='\t',header=None,low_memory=False)
            Zscore.columns = np.array(tuple(self.Zscore_header))
            
            Weight = pd.read_csv(StringIO(Weight_out.decode('utf-8')),sep='\t',header=None,low_memory=False)
            Weight.columns = np.array(tuple(self.Weight_header))
            
            Zscore['snpID'] = (Zscore['CHROM'].astype('str')+':'+Zscore['POS'].astype('str')
                               +':'+Zscore.REF+':'+Zscore.ALT)
            
            Weight['snpID'] = (Weight['CHROM'].astype('str')+':'+Weight['POS'].astype('str')
                               +':'+Weight.REF+':'+Weight.ALT)

            snp_overlap = np.intersect1d(np.array(Weight.snpID),np.array(Zscore.snpID))
            # Merge Z-score and Weight File
            ZW = (Weight >> mask(Weight.snpID.isin(snp_overlap))).merge((Zscore
                                                                         >>mask(Zscore.snpID.isin(snp_overlap))
                                                                         >>select(Zscore.snpID,Zscore.Zscore)),
                                                                        left_on='snpID',
                                                                        right_on='snpID')
            ZW = ZW.drop_duplicates(['snpID'],keep='first')
            ZW = ZW.sort_values(by='POS').reset_index(drop=True)

            ### Read in Reference Covariance Matrix File
            covar_process = subprocess.Popen(["tabix"+" "+self.Covar+" "+str(self.chr_num)+":"+str(min(ZW.POS))+"-"+str(max(ZW.POS))],
                                             shell=True,
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
            covar_out=covar_process.communicate()[0]

            if len(covar_out) == 0:
                print("No Reference Covariance Information for Gene:"+self.TargetID[num])
            else:
                MCOV = pd.read_csv(StringIO(covar_out.decode('utf-8')),sep='\t',header=None,low_memory=False)
                MCOV.columns = ['CHROM','POS','ID','REF','ALT','COV']
                
                MCOV['snpID'] = (MCOV['CHROM'].astype('str')+':'+MCOV['POS'].astype('str')
                                 +':'+MCOV.REF+':'+MCOV.ALT)
                MCOV = MCOV.drop_duplicates(['snpID'],keep='first')
                
                MCOV['COV'] = MCOV['COV'].apply(lambda x:np.array(x.split(',')).astype('float'))
                MCOV = MCOV.sort_values(by='POS').reset_index(drop=True)
                indicates = (MCOV >> mask(MCOV.snpID.isin(ZW.snpID))).index
                
                V_temp = np.zeros((len(indicates),len(indicates)))
                
                for i in range(len(indicates)):
                    cov_temp = MCOV.COV.loc[indicates[i]]
                    N = len(cov_temp)
                    
                    for j in range(i,len(indicates)):
                        if indicates[j] - indicates[i] < N:
                            V_temp[i,j] = cov_temp[indicates[j]-indicates[i]]
                        else:
                            V_temp[i,j] = 0
                V=V_temp+V_temp.T-np.diag(V_temp.diagonal())
                
                ZW = ZW >> mask(ZW.snpID.isin(MCOV.snpID))
                ### Calculate Burden Z score
                burden_Z = mat(ZW.Zscore)*mat(ZW.ES).T/sqrt(mat(ZW.ES)*V*mat(ZW.ES).T)
                
                result = pd.DataFrame()
                result['CHROM'] = np.array(self.chr_num).ravel()
                result['GeneStart'] = np.array(Gene_temp.GeneStart).ravel()
                result['GeneEnd'] = np.array(Gene_temp.GeneEnd).ravel()
                result['TargetID'] = np.array(self.TargetID[num]).ravel()
                result['GeneName'] = np.array(Gene_temp.GeneName).ravel()
                result['ZSCORE'] = burden_Z
                ### p-value for chi-square test
                result['PVALUE'] = 1-chi2.cdf(burden_Z**2,1)
                
                result.to_csv(self.out_prefix+'/CHR'+str(self.chr_num)+'_association_study.txt',
                              sep='\t',index=None,header=None,mode='a')
                
    def main(self):
        pool = multiprocessing.Pool(self.thread) 

        pool.map(self.thread_process,[num for num in range(len(self.TargetID))])
        pool.close()
        pool.join()
                
                










