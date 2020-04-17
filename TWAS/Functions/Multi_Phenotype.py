#!/usr/bin/env python
###########################################################################
### import package need
import pandas as pd
import numpy as np
from numpy import *
from dfply import *
from scipy import stats
import multiprocessing

# For OLS and Logistics regression
import statsmodels.api as sm
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

###############################################################################################
### Input:
### 1. method: Link Function (OLS Only).
### 2. pheno: An array of Phenotype Names
### 3. covar: An array of Covariance Names
### 4. PED: PED File
### 5. Gene_Exp: Gene Annotation and Expression Level File
### 6. TargetID
### 7. out_prefix: Output Dir
### 8. thread: Number of Thread

### Output:
### 1. CHROM, GeneStart, GeneEnd , GeneName , TargetID: Gene Annotation columns
### 2. R2: Regression R-square
### 3. F_STAT: Value of F statistics of Regression Model
### 4. F_PVALUE: P-value of F Test of Regression Model
### 5. N: Sample Size

###############################################################################################
# For Multiple Phenotype
##############################################################################################
class Multi_Phenotype(object):
    def __init__(self, method, pheno, covar, PED, Gene_Exp, 
                 TargetID, out_prefix, thread):
        self.method = method
        self.pheno = pheno
        self.covar = covar
        self.PED = PED
        self.Gene_Exp = Gene_Exp
        self.TargetID = TargetID
        self.out_prefix = out_prefix
        self.thread = thread
    
    def regression_multi(self, X, Y, TargetID):
        lm = sm.OLS(Y,X).fit()
        
        result = pd.DataFrame()
        result['TargetID'] = np.array(TargetID).ravel()
        # Regression R-square
        result['R2'] = np.array(lm.rsquared).ravel()
        # F Statistics
        result['F_STAT'] = np.array(lm.fvalue).ravel()
        # F Test P-value
        result['F_PVALUE'] = np.array(lm.f_pvalue).ravel()
        # Sample Size for Regression
        result['N'] = np.array(len(X))
        
        return result
    
    def thread_process(self, num):
        Target_temp = self.Target >> select(self.Target[self.TargetID[num]],self.Target[self.pheno])
        Target_temp = pd.DataFrame(Target_temp.dropna(axis=0,how='any'),dtype='float')
        
        X = Target_temp[self.pheno]
        Y = Target_temp[self.TargetID[num]]
        
        lm = self.regression_multi(X,Y,self.TargetID[num])
        
        Gene_annot = self.Gene_Exp >> mask(self.Gene_Exp.TargetID==self.TargetID[num]) >> select(self.Gene_Exp.columns[0:5])
        out = Gene_annot.merge(lm,left_on='TargetID',right_on='TargetID',how='outer')
        
        out.to_csv(self.out_prefix+"/association_study_Multi_"+self.method+".txt",
                   sep='\t',header=None,index=None,mode='a')
    
    def main(self):
        Gene_temp = (self.Gene_Exp >> select(self.Gene_Exp.TargetID, 
                                             self.Gene_Exp[self.Gene_Exp.columns[5:]])).T
        Gene_temp.columns = Gene_temp.loc['TargetID']
        Gene_temp = Gene_temp.drop(['TargetID'])
        Gene_temp['IND_ID'] = Gene_temp.index
        Gene_temp = Gene_temp.reset_index(drop=True)

        ### Residuals Calculation
        res = pd.DataFrame()
        res['IND_ID'] = np.array(self.PED.IND_ID).ravel()
        
        for i in range(len(self.pheno)):
            res[self.pheno[i]] = np.array(sm.OLS(self.PED[self.pheno[i]],
                                                 sm.add_constant(self.PED[self.covar])).fit().resid).ravel()
        ### Initializing Header of Output Files
        pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
                              'R2','F_STAT','F_PVALUE','N']).to_csv(self.out_prefix+"/association_study_Multi_"+self.method+".txt",
                                                                    sep='\t',header=True,index=None,mode='a')
        
        Target = res.merge(Gene_temp,left_on='IND_ID',right_on='IND_ID',how='outer')
        Target = pd.DataFrame((Target >> drop(Target.IND_ID)),dtype='float')

        self.Target = Target

        pool = multiprocessing.Pool(self.thread)
        pool.map(self.thread_process,[num for num in range(len(self.TargetID))])
        pool.close()
        pool.join()







