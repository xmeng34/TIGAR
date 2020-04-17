#!/usr/bin/env python
###############################################################################
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

###############################################################################
### Input:
### 1. method: Link Function. OLS or Logit.
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
### 3. BETA: Regression Coefficient of Gene
### 4. BETA_SE: Standard Error of BETA
### 5. T_STAT: Value of T statistics of Gene
### 6. PVALUE: P-value of T Test of Gene
### 7. N: Sample Size

################################################################################
# For Single Phenotype
################################################################################
class Single_Phenotype(object):
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
    
    def regression_single(self,method,X,Y,TargetID):
        ### Add Intercept Columns for Design Matrix
        newX = sm.add_constant(X)
        
        result = pd.DataFrame()
        result['TargetID'] = np.array(TargetID).ravel()
        
        ### Regression
        if method=='OLS':
            lm = sm.OLS(Y, newX).fit()
            # Regression R-square
            result['R2'] = np.array(lm.rsquared)
        
        elif method=='Logit':
            lm = sm.Logit(Y-1, newX).fit(maxiter=1000)
            # Regression R-square
            result['R2'] = np.array(lm.prsquared)

        # Regression coefficient for Gene
        result['BETA'] = pd.DataFrame(lm.params).loc[TargetID]
        # Standard Error for Gene
        result['BETA_SE'] = pd.DataFrame(lm.bse).loc[TargetID]
        # T Statistics for Gene
        result['T_STAT'] = pd.DataFrame(lm.tvalues).loc[TargetID]
        # T Test p-value for Gene
        result['PVALUE'] = pd.DataFrame(lm.pvalues).loc[TargetID]
        # Sample Size for Regression
        result['N'] = np.array(len(X))

        return result

    def thread_process(self, num):
        Target_temp = (self.Target 
                       >> select(self.Target[self.pheno],self.Target[self.covar],
                                 self.Target[self.TargetID[num]])).dropna(axis=0,how='any')
        X = (self.Target >> select(self.Target[self.covar],self.Target[self.TargetID[num]]))
        Y = self.Target[self.pheno]
    
        lm = self.regression_single(self.method,X,Y,self.TargetID[num])
        
        Gene_annot = self.Gene_Exp >> mask(self.Gene_Exp.TargetID==self.TargetID[num]) >> select(self.Gene_Exp.columns[0:5])
    
        out = Gene_annot.merge(lm,left_on='TargetID',right_on='TargetID',how='outer')

        out.to_csv(self.out_prefix+"/association_study_Single_"+self.method+".txt",sep='\t',header=None,index=None,mode='a')

    def main(self):
        Gene_temp = (self.Gene_Exp >> select(self.Gene_Exp.TargetID, 
                                             self.Gene_Exp[self.Gene_Exp.columns[5:]])).T

        Gene_temp.columns = Gene_temp.loc['TargetID']
        Gene_temp = Gene_temp.drop(['TargetID'])
        Gene_temp['IND_ID'] = Gene_temp.index
        Gene_temp = Gene_temp.reset_index(drop=True)
        
        Target = self.PED.merge(Gene_temp, left_on='IND_ID', right_on='IND_ID', how='outer')
        Target = pd.DataFrame((Target >> drop(Target.IND_ID)),dtype='float')

        self.Target = Target
        
        ### Initializing Header of Output Files
        pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
                              'R2','BETA','BETA_SE','T_STAT','PVALUE','N']).to_csv(self.out_prefix+"/association_study_Single_"+self.method+".txt",
                                                                                   sep='\t',header=True,index=None,mode='a')
        
        pool = multiprocessing.Pool(self.thread)
        pool.map(self.thread_process,[num for num in range(len(self.TargetID))])
        pool.close()
        pool.join()



