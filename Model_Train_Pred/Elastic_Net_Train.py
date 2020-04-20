#!/usr/bin/env python

##############################################################################################
# Import Packages Needed
import argparse
import warnings
from time import *
import subprocess
from subprocess import *
import io
from io import StringIO
from io import *
import pandas as pd 
import numpy as np 
from numpy import *
from dfply import *
import multiprocessing

### Import Grid Search for Model Selection
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold

### For Elastic Net Regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import r2_score

### For OLS Regression In Cross Validation
import statsmodels.api as sm
from scipy import stats

import sys
sys.path.append("./Model_Train_Pred/Functions")
from CHR_Reform import CHR_Reform

###############################################################################################
# Hide Warnings from Python
warnings.filterwarnings("ignore")
###############################################################################################
### Time Calculation
start_time=time()

###############################################################################################
# Parameter Selection & Elastic Net Regression
#############################################################################################
### Elastic Net
### Input:
### 1.trainX: Values for Independent Variables (SNPs Genotype)
### 2.trainY: Values for Response Variable (Gene Expression Level)
### 3.testX: Values for Independent Variables for Prediction
### 4.testY: Values for Response Variables for Prediction
### 5.k: K-Fold Cross Validation
### 6.Alpha: Ratio for L1 and L2 Penalty in Elastic Net Regression
###          Alpha=0: Lasso Regression
###          Alpha=1: Ridge Regression
###          0 < Alpha < 1: Elastic Net Regression

### Return:
### 1.Regression Coefficients
### 2.Training R-square
### 3.Parameter Selected by Cross Validation
### 4.Mean Cross Validation Score

### Using Grid Search and Cross Validation to Find the Best Lambda (Penalty Term)

def elastic_net(trainX,trainY,testX,testY,k,Alpha):
    clf = GridSearchCV(ElasticNet(l1_ratio=Alpha,fit_intercept=False),
                       [{'alpha':np.arange(0,1.01,0.01)}],cv=k).fit(trainX,trainY)
    
    reg = ElasticNet(l1_ratio=Alpha,alpha=clf.best_params_['alpha']).fit(trainX,trainY)

    lm = sm.OLS(testY,sm.add_constant(reg.predict(testX))).fit()

    return reg.coef_,lm.rsquared,lm.f_pvalue,clf.best_params_['alpha'],clf.best_score_


##############################################################################################
### Create Calculation Matrix
def Info_Gene(Chr, Exp):
    sampleID = np.intersect1d(Chr.columns, Exp.columns[5:])
    
    CHR_target = Chr >> select(Chr.snpID, Chr[sampleID])
    Exp_target = Exp >> select(Exp.TargetID, Exp[sampleID])
    
    info_temp = pd.concat([CHR_target,Exp_target],sort=False).T
    
    snpID = np.array(info_temp.loc['snpID'].dropna())
    
    TargetID = Exp.TargetID
    
    info_temp.columns = np.hstack((snpID,TargetID))
    
    info = pd.DataFrame(info_temp.drop(['snpID','TargetID']),dtype=np.float)
    
    return info,snpID, sampleID

###############################################################################################
# Model Training
###############################################################################################
### Variables Need
parser = argparse.ArgumentParser(description='Inputs for this Script.')

### Gene Annotation and Expression Level File
parser.add_argument('--Gene_Exp_path',type=str,default = None)

### Training SampleIDs
parser.add_argument('--sampleID', type=str, default = None)

### Specified Chromosome Number
parser.add_argument('--chr_num',type=int,default = None)

### Genotype File
parser.add_argument('--genofile', type=str, default=None)
parser.add_argument('--genofile_header', type=str, default=None)

### Specified Input Genotype File Type (vcf or dosages)
parser.add_argument('--geno',type=str, choices=['vcf', 'dosages'], default = None)

### GT or DS
parser.add_argument('--Format',choices=['GT', 'DS'],type=str,default=None)

### Threshold for Data Selection
### Minor Allele Frequency (Range from 0-1)
parser.add_argument('--maf',type=float,default=None)

### p-value for Hardy Weinberg Equilibrium Exact Test
parser.add_argument('--hwe',type=float,default=None)

### Window Size around Gene Boundary
parser.add_argument('--window',type=int,default=None)

### Parameter Selection for Model Training
### K-Fold Cross Validation
parser.add_argument('--cv',type=int,default=None)
### Ratio of L1 and L2
parser.add_argument('--alpha',type=float,default=None)

### Number of Thread
parser.add_argument('--thread',type=int,default = None)

### Output Dir
parser.add_argument('--out_prefix',type=str,default=None)

args = parser.parse_args()

##############################################################################################
### Checking Input Commands
print("Gene Annotation and Expression Level File:"+args.Gene_Exp_path)
print("Training SampleIDs:"+args.sampleID)
print("Chromosome Number:"+str(args.chr_num))
print("Training Genotype File:"+args.genofile)
print(args.genofile_header)

if args.geno=='vcf':
    print("Using "+args.Format+" Format for Training")
elif args.geno=='dosages':
    print("Using DS Format for Training.")

print("Threshold for MAF:"+str(args.maf))
print("Threshold for p-value of HW Test:"+str(args.hwe))
print("window="+str(args.window))
print("Using "+str(args.cv)+"-fold for Cross-Validation")
print("L1 & L2 ratio in Elastic Net Regression:"+str(args.alpha))
print("Number of Thread:"+str(args.thread))
print("Output Dir:"+args.out_prefix)

###############################################################################################
### Training Processing

### Read in Gene Annotation and Expression Level File (.txt file)
### First Five Columns should be Fixed:
### 1.Chromosome Number (CHROM)
### 2.Gene Starting Position (GeneStart)
### 3.Gene Ending Position (GeneEnd)
### 4.TargetID (i.e. GeneID, treated as unique annotation for each gene)
### 5.Gene Name (GeneName)

Gene_Exp = pd.read_csv(args.Gene_Exp_path,sep='\t',low_memory=False)
Gene_header = np.array(Gene_Exp.columns)
Gene_header[0:5] = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName']
Gene_Exp.columns = Gene_header
if len(Gene_Exp)==0:
    raise SystemExit("No Gene and Expression Level Data for this Chromosome.")

### SampleIDs
sampleID_temp = pd.read_csv(args.sampleID,sep='\t',header=None)
sampleID_temp = np.array(sampleID_temp).ravel()

### Initializing Header of Output Files
pd.DataFrame(columns=['CHROM','POS','REF','ALT','TargetID', 
                      'ID', 
                      'ES','MAF','p_HWE']).to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_elastic_net_training_weight.txt',
                                                  header=True,index=None,sep='\t',mode='a')

pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
                      'snp_size','effect_snp_size','sample_size',
                      '5-fold-CV-R2','TrainPVALUE','Train-R2',
                      'k-fold','alpha','Lambda','cvm']).to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_elastic_net_training_info.txt',
                                                               header=True,index=None,sep='\t',mode='a')

### Read in Header for Genotype File
genofile_header=pd.read_csv(args.genofile_header,sep='\t').rename(columns={'#CHROM':'CHROM'})
genofile_header=np.array(tuple(genofile_header))
#############################################################################################
### Thread Process
def thread_process(num):
    Exp_temp = pd.DataFrame(Gene_Exp.loc[num]).T
    TargetID = np.array(Exp_temp.TargetID)[0]
    
    ### Gene Boundary
    start = max(int(Exp_temp.GeneStart)-args.window,0)
    end = max(int(Exp_temp.GeneEnd)+args.window,0)
    
    ### Select Corresponding Genotype File by Tabix
    chr_process = subprocess.Popen(["tabix"+" "+args.genofile+" "+str(args.chr_num)+":"+str(start)+"-"+str(end)],
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
    out = chr_process.communicate()[0]
    
    if len(out) == 0:
        print("No Corresponding Genotype Data for this Gene:"+TargetID)
    
    else:
        ### Decode subprocess Output with 'utf-8'
        Chr_temp = pd.read_csv(StringIO(out.decode('utf-8')),sep='\t',header=None,low_memory=False)
        Chr_temp.columns = genofile_header
        Chr_temp = Chr_temp.reset_index(drop=True)
        
        ### Reform Genotype Data
        Chrom =  CHR_Reform('TRAIN', Chr_temp, args.geno, sampleID_temp, args.Format, p_hwe=args.hwe, maf=args.maf).main()
        
        Information, SNPs, sampleID = Info_Gene(Chrom, Exp_temp)

        sample_size=len(sampleID)

        ### 5-Folds Cross Validation before Training Start
        ### Separate sampleIDs for 5-Folds Cross Validation
        CV_trainID = []
        CV_testID = []

        kf=KFold(n_splits=5)
        for train_index,test_index in kf.split(sampleID):
            CV_trainID.append(np.array(','.join(sampleID[train_index])))
            CV_testID.append(np.array(','.join(sampleID[test_index])))

        CV_trainID = pd.DataFrame(CV_trainID).apply(lambda x:x.str.split(","))
        CV_testID = pd.DataFrame(CV_testID).apply(lambda x:x.str.split(","))

        ### Evaluate whether Elastic Net Model Valid for this Gene
        k_fold_R2 = []
        for i in range(5):
            Info_CV_train=Information.loc[np.array(CV_trainID.loc[i][0])].dropna()
            Info_CV_test=Information.loc[np.array(CV_testID.loc[i][0])].dropna()

            k_fold_R2.append(elastic_net(Info_CV_train[SNPs],Info_CV_train[TargetID],
                                         Info_CV_test[SNPs],Info_CV_test[TargetID],
                                         args.cv,args.alpha)[1])

        flag = sum(k_fold_R2)/5
        
        if flag < 0.01:
            print("Elastic Net model is not Valid for Gene:"+TargetID)
            print("5-fold-CV-R2="+str(flag))
            
        else:
            print("Running Elastic Net Training from Gene:"+TargetID)
            print("5-fold-CV-R2="+str(flag))
            ### Store Training Parameters
            Beta_temp = pd.DataFrame()
            Beta_temp['SNPs'] = SNPs
            Beta_temp['TargetID'] = TargetID
            Beta_temp['beta'], R2, Pvalue, Lambda, cvm = elastic_net(Information[SNPs],Information[TargetID], 
                                                                     Information[SNPs], Information[TargetID],
                                                                     args.cv,args.alpha)
            Beta=(Chrom[['CHROM','POS','ID','REF','ALT','snpID',
                         'p_HWE','MAF']].merge(Beta_temp,left_on='snpID',
                                               right_on='SNPs',how='outer'))

            Beta = Beta >> drop(Beta[['SNPs', 'snpID']])
            
            ### Only Keep snps with beta!=0
            Beta = (Beta
                    >> mask(Beta.beta != 0) 
                    >> select(Beta.CHROM,Beta.POS,Beta.REF,Beta.ALT,Beta.TargetID,
                              Beta.ID,Beta.beta,Beta.MAF,Beta.p_HWE))
            Beta = Beta.rename(columns = {'beta':'ES'})
            
            ### Training Weight
            Beta.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_elastic_net_training_weight.txt',
                        header=None, index=None,sep='\t',mode='a')

            ### Store Training Information
            ### Store Result from Elastic Net
            Info_Train = pd.DataFrame()

            Info_Train['TargetID'] = np.array(TargetID).ravel()
            Info_Train['snp_size'] = len(SNPs)
            Info_Train['effect_snp_size'] = len(Beta.ES)
            Info_Train['sample_size'] = sample_size
            Info_Train['5-fold-CV-R2'] = np.array(flag).ravel()
            Info_Train['TrainPVALUE'] = np.array(Pvalue).ravel()
            
            if len(Beta.ES) == 0:
                Info_Train['Train-R2'] = 0
            else:
                Info_Train['Train-R2'] = R2

            Info_Train['k_fold']=args.cv
            Info_Train['alpha']=args.alpha
            Info_Train['Lambda']=Lambda
            Info_Train['cvm']=cvm
            
            Info = (Exp_temp >> select(Exp_temp[['CHROM','GeneStart','GeneEnd',
                                                 'GeneName','TargetID']])).merge(Info_Train,
                                                                                 left_on='TargetID',
                                                                                 right_on='TargetID',
                                                                                 how='outer')
            ### Training Information
            Info.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_elastic_net_training_info.txt',
                        header=None, index=None,sep='\t',mode='a')

####################################################################################################
### Start Thread
if (args.thread < int(len(Gene_Exp)/100) | args.thread > len(Gene_Exp)):
    args.thread = (int(len(Gene_Exp)/100)+1)*100

pool = multiprocessing.Pool(args.thread)

pool.map(thread_process,[num for num in range(len(Gene_Exp))])

pool.close()
pool.join()

#####################################################################################################
### Time Calculation
end_time=time()
print("Running Time:"+str(round((end_time-start_time)/60,2))+" minutes.")






