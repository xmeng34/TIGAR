#!/usr/bin/env python

###############################################################################################
# Import Packages Needed
import argparse
from time import *
import subprocess
from subprocess import *
import shlex
import io
from io import StringIO
import pandas as pd
import numpy as np
from numpy import *
from dfply import *

from sklearn.model_selection import KFold
import statsmodels.api as sm
import scipy.stats as stats

import multiprocessing

import sys
sys.path.append("./Model_Train_Pred/Functions")
from CHR_Reform import CHR_Reform

#################################################################################################
# Hide Warnings from Python
warnings.filterwarnings("ignore")
#################################################################################################
### Time Calculation
start_time=time()

#################################################################################################
# Dealing with DPR Input Files
#################################################################################################
### Input:
### 1. data : DataFrame for extracting DPR input files.
### 2. sampleID : sampleIDs contained in bimbam, pheno and SNP annotation file
### 3. TargetID : Specifying Gene for running DPR

### Output:
### 1. bimbam file
### 2. phenotype file
### 3. SNP annotation file

def DPR_Input_Generator(data, sampleID, TargetID):
    DPR_Info = data >> select(data.snpID, data.POS, data.CHROM, data.REF, data.ALT, 
                              data[sampleID], data.TargetID)
    
    bimbam = (DPR_Info >> select(DPR_Info.snpID,
                                 DPR_Info.REF,
                                 DPR_Info.ALT,
                                 DPR_Info[sampleID])).dropna(axis=0,how='any')
    
    pheno = (DPR_Info 
             >> mask(DPR_Info.TargetID==TargetID) 
             >> drop(DPR_Info.TargetID)).dropna(axis=1,how='any')
    
    SNP_annot = (DPR_Info 
                 >> select(DPR_Info.snpID,DPR_Info.POS,DPR_Info.CHROM)).dropna(axis=0,how='any')
    
    return bimbam, pheno, SNP_annot

#################################################################################################
# Model Training
#################################################################################################
### Variables Need
parser = argparse.ArgumentParser(description='Inputs for this Scripts')

### for Gene Annotation and Expression Level File
parser.add_argument('--Gene_Exp_path',type=str,default = None)

### Training SampleIDs
parser.add_argument('--sampleID',type=str,default = None)

### Specified Chromosome Number
parser.add_argument('--chr_num',type=int,default = None)

### Genotype File
parser.add_argument('--genofile',type=str,default = None)
parser.add_argument('--genofile_header',type=str,default = None)

### Specified Input Genotype File Type (vcf or dosages)
parser.add_argument('--geno',choices=['vcf', 'dosages'],type=str,default = None)

### GT or DS
parser.add_argument('--Format',choices=['GT', 'DS'],type=str,default=None)

### Threshold for Data Selection
### Minor Allele Frequency (Range from 0-1)
parser.add_argument('--maf',type=float,default=None)

### p-value for Hardy Weinberg Equilibrium Exact Test
parser.add_argument('--hwe',type=float,default=None)

### Window Size around Gene Boundary
parser.add_argument('--window',type=int,default=None)

### Model to Run DPR
parser.add_argument('--dpr',choices=[1,2,3],type=int,default=None)

### Define Effect-Size
### fixed : ES=beta
### additive : ES=beta+b
parser.add_argument('--ES',choices=['fixed', 'additive'],type=str,default=None)

### Number of Thread
parser.add_argument('--thread',type=int,default = None)

### Output Dir
parser.add_argument('--out_prefix',type=str,default=None)

args = parser.parse_args()

#######################################################################################################################
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
print("Threshold for p-value of HW test:"+str(args.hwe))
print("window="+str(args.window))
print("Running DPR Model with dpr="+str(args.dpr))
print("Using Effect-Size:"+args.ES)
print("Number of Thread:"+str(args.thread))
print("Output Dir:"+args.out_prefix)
#######################################################################################################################
# Training Processing

### Read in Gene Annotation and Expression Level File (.txt file)
### First Five Columns should be Fixed:
### 1.Chromosome Number (CHROM)
### 2.Gene Starting Position (GeneStart)
### 3.Gene Ending Position (GeneEnd)
### 4.TargetID (i.e.GeneID, treated as unique annotation for each gene)
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
                      'n_miss','b','beta','gamma',
                      'ES','MAF','p_HWE']).to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_DPR_training_weight.txt',
                                                  sep='\t',header=True,index=None,mode='a')

pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
                      'snp_size','effect_snp_size','sample_size',
                      '5-fold-CV-R2','TrainPVALUE','Train-R2']).to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_DPR_training_info.txt',
                                                                       sep='\t',header=True,index=None,mode='a')

### Read in Header for Genotype File
genofile_header=pd.read_csv(args.genofile_header,sep='\t').rename(columns={'#CHROM':'CHROM'})
genofile_header=np.array(tuple(genofile_header))

#############################################################################################
### Thread Process
def thread_process(num):
    Exp_temp = pd.DataFrame(Gene_Exp.loc[num]).T
    TargetID = np.array(Exp_temp.TargetID)[0]
    
    ### Gene Boundary
    start=max(int(Exp_temp.GeneStart)-args.window,0)
    end=max(int(Exp_temp.GeneEnd)+args.window,0)

    ### Select Corresponding Genotype File by Tabix
    chr_process=subprocess.Popen(["tabix"+" "+args.genofile+" "+str(args.chr_num)+":"+str(start)+"-"+str(end)],
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
    out=chr_process.communicate()[0]

    if len(out)==0:
        print("No Corresponding Genotype Data for this Gene:"+TargetID)

    else:
        print("Preparing DPR Input for Gene:"+TargetID)
        ### Decode subprocess Output with 'utf-8'
        Chr_temp=pd.read_csv(StringIO(out.decode('utf-8')),sep='\t',header=None,low_memory=False)
        Chr_temp.columns = genofile_header
        Chr_temp=Chr_temp.reset_index(drop=True)
        
        ### Reform Genotype Data
        Chrom =  CHR_Reform('TRAIN', Chr_temp, args.geno, sampleID_temp, args.Format, p_hwe=args.hwe, maf=args.maf).main()
        
        ### Overlap sampleIDs in Genotype and Gene Expression Data
        sampleID = np.intersect1d(Chrom.columns, Exp_temp.columns[5:])
        sample_size=len(sampleID)

        Info = pd.concat([Chrom >> select(Chrom[['CHROM', 'POS', 'REF', 'ALT', 'snpID']], Chrom[sampleID]),
                          Exp_temp >> select(Exp_temp.TargetID, Exp_temp[sampleID])],
                         sort=False)
        
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

        ### Evaluate whether DPR Model Valid for this Gene
        CV_file_dir=args.out_prefix+"/DPR_input/CV"

        k_fold_R2=[]
        k_fold_R2.append(0)
        for i in range(5):
            trainID = np.array(CV_trainID.loc[i][0])
            testID = np.array(CV_testID.loc[i][0])

            ### Create DPR Input File for Train Set
            bimbam_train, pheno_train, SNP_annot_train = DPR_Input_Generator(Info, trainID, TargetID)

            ### bimbam file
            bimbam_train.to_csv(CV_file_dir+'/bimbam/'+TargetID+'_CV'+str(i+1)+'_bimbam.txt',
                                header=False,index=None,sep='\t')
            ### phenotype file
            pheno_train.T.to_csv(CV_file_dir+'/pheno/'+TargetID+'_CV'+str(i+1)+'_pheno.txt',
                                 header=False,index=None,sep='\t')
            ### SNP annotation file
            SNP_annot_train.to_csv(CV_file_dir+'/SNP_annot/'+TargetID+'_CV'+str(i+1)+'_snp_annot.txt',
                                   header=False,index=None,sep='\t')
            ### call DPR
            TargetID_CV = TargetID+'_CV'+str(i+1)
            stop_CV=0
            try:
                subprocess.check_call(shlex.split('./Model_Train_Pred/Functions/call_DPR.sh'+' '+CV_file_dir+' '+str(args.dpr)+' '+TargetID_CV))
            except subprocess.CalledProcessError as err:
                stop_CV=1
                print("DPR failed in CV"+str(i+1)+" for TargetID:"+TargetID)
                
            if stop_CV==1:
                continue
            else:
                ### Read in DPR Results for Cross Validation
                result_CV = pd.read_csv(CV_file_dir+'/output/DPR_'+TargetID_CV+'.param.txt',sep='\t').rename(columns={'rs':'snpID'})

                ### Overall Effect Size
                if args.ES=='fixed':
                    result_CV['ES'] = result_CV['beta']
                elif args.ES=='additive':
                    result_CV['ES'] = result_CV['b']+result_CV['beta']
                else:
                    raise SystemExit("Effect-size can not identify.")

                result_CV = result_CV >> select(result_CV.snpID,result_CV.ES)
            
                ### Calculate Predicted Value for Test Set
                bimbam_test, pheno_test = DPR_Input_Generator(Info, testID, TargetID)[0:2]
                pheno_test = pheno_test.reset_index(drop=True)
                pheno_test = pd.DataFrame(pheno_test,dtype=np.float) 

                overall = bimbam_test.merge(result_CV,left_on='snpID',right_on='snpID',how='outer')
                overall['ES'].fillna(0,inplace=True)

                Ypred=np.array(mat(pd.DataFrame(overall[testID],dtype=np.float)).T*mat(overall.ES).reshape((len(overall.snpID),1))).ravel()

                lm = sm.OLS(np.array(pheno_test.loc[0]),sm.add_constant(Ypred)).fit()
                k_fold_R2.append(lm.rsquared)

        flag = sum(k_fold_R2)/5
        
        if flag < 0.01:
            print("DPR model is not Valid for Gene:"+TargetID)
            print("5-fold-CV-R2="+str(flag))
        
        else:
            print("Running DPR Training for Gene:"+TargetID)
            print("5-fold-CV-R2="+str(flag))

            ### Create DPR Input Files for Training
            file_dir=args.out_prefix+'/DPR_input'
            print("Running model training for Gene:"+TargetID)
            bimbam, pheno, SNP_annot = DPR_Input_Generator(Info, sampleID, TargetID)
            bimbam.to_csv(file_dir+'/bimbam/'+TargetID+'_bimbam.txt',
                          header=False,index=None,sep='\t')
            pheno.T.to_csv(file_dir+'/pheno/'+TargetID+'_pheno.txt',
                           header=False,index=None,sep='\t')
            SNP_annot.to_csv(file_dir+'/SNP_annot/'+TargetID+'_snp_annot.txt',
                             header=False,index=None,sep='\t')
            
            stop_DPR=0
            try:
                subprocess.check_call(shlex.split('./Model_Train_Pred/Functions/call_DPR.sh'+' '+file_dir+' '+str(args.dpr)+' '+TargetID))
            except subprocess.CalledProcessError as err:
                stop_DPR=1
                print("DPR failed for TargetID:"+TargetID)
            
            if stop_DPR==0:
                Info_Train=pd.DataFrame()
                Info_Train['TargetID']=np.array(TargetID).ravel()
                
                result=pd.read_csv(file_dir+'/output/DPR_'+TargetID+'.param.txt',sep='\t')
                result['TargetID']=TargetID

                if args.ES=='fixed':
                    result['ES']=result['beta']
                elif args.ES=='additive':
                    result['ES']=result['beta']+result['b']

                Info_Train['snp_size']=np.array(len(result.ES)).ravel()
                ### Only Keep snps with ES!=0
                result = result >> mask(result.ES!=0)
                Info_Train['effect_snp_size']=np.array(len(result.ES)).ravel()

                ### Output DPR Training Result to a Single File
                result=result.rename(columns={'chr':'CHROM','rs':'snpID','ps':'POS'})

                result=result >> select(result.CHROM,result.snpID,result.POS,result.TargetID,
                                        result.n_miss,result.b,result.beta,result.ES,result.gamma)

                ### Store Filter Information
                Filter = Chrom >> select(Chrom.snpID,Chrom.p_HWE,Chrom.MAF)
        
                param = result.merge((Filter
                                      >> mask(Filter.snpID.isin(np.array(result.snpID)))),
                                     left_on='snpID',right_on='snpID',how='outer')
                param['REF'] = param['snpID'].apply(lambda x:x.split(":")[2])
                param['ALT'] = param['snpID'].apply(lambda x:x.split(":")[3])

                param = param >> select(param.CHROM,param.POS,param.REF,param.ALT,param.TargetID,
                                        param.n_miss,param.b,param.beta,param.gamma,param.ES,param.MAF,param.p_HWE)
                param['CHROM'] = param['CHROM'].astype('int')
                param['POS'] = param['POS'].astype('int')

                ### Training Weight
                param.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_DPR_training_weight.txt',
                             sep='\t',header=None,index=None,mode='a')

                result = result >> select(result.snpID,result.ES)

                ### for R2 calculation
                bimbam = bimbam >> drop(bimbam.REF,bimbam.ALT)

                ### Read in Phenotype File
                pheno = pd.DataFrame(pheno,dtype=np.float).reset_index(drop=True)

                ID = np.intersect1d(np.array(result.snpID),np.array(bimbam.snpID))

                pred_temp = (result
                             >> mask(result.snpID.isin(ID))).merge((bimbam
                                                                    >> mask(bimbam.snpID.isin(ID))),
                                                                   left_on='snpID',right_on='snpID',
                                                                   how='outer')
                pred = pred_temp.T
                pred.columns = pred.loc['snpID']
                pred = (pred.drop(['snpID','ES'])).reset_index(drop=True)
                pred = pd.DataFrame(pred,dtype='float')
                pheno_pred = np.array(mat(pred)*mat(pred_temp.ES).T).ravel()

                lm_final = sm.OLS(np.array(pheno_pred),sm.add_constant(np.array(pheno.loc[0]))).fit()
                Info_Train['sample_size'] = sample_size
                Info_Train['5-fold-CV-R2'] = np.array(flag).ravel()
                Info_Train['TrainPVALUE'] = np.array(lm_final.f_pvalue).ravel()
                Info_Train['Train-R2'] = np.array(lm_final.rsquared).ravel()

                Info = (Exp_temp
                        >> select(Exp_temp[['CHROM','GeneStart','GeneEnd',
                                            'GeneName','TargetID']])).merge(Info_Train,
                                                                            left_on='TargetID',
                                                                            right_on='TargetID',
                                                                            how='outer')
                ### Training Information
                Info.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_DPR_training_info.txt',
                            sep='\t',header=None,index=None,mode='a')

################################################################################################################
### Start Thread
if (args.thread < int(len(Gene_Exp)/100) | args.thread > len(Gene_Exp)):
    args.thread = (int(len(Gene_Exp)/100)+1)*100

pool = multiprocessing.Pool(args.thread)

pool.map(thread_process,[num for num in range(len(Gene_Exp))])

pool.close()
pool.join()
#################################################################################################################
### Time Calculation
end_time=time()
print("Running Time:"+str(round((end_time-start_time)/60,2))+" minutes.")














