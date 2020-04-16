#!/usr/bin/env python

##########################################################################################
# Import Packages Needed
import argparse
from time import *
import os
import subprocess
from subprocess import *
import io
from io import StringIO
import pandas as pd
import numpy as np
from numpy import *
from dfply import *
import multiprocessing

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
# Model Prediction
###############################################################################################
### Variables Need
parser = argparse.ArgumentParser(description='Inputs for this Script.')

### Specified Training Model
parser.add_argument('--model',choices=['elastic_net', 'DPR'],type=str,default = None)

### Specified Chromosome Number
parser.add_argument('--chr_num',type=int,default = None)

### Training Weight File
parser.add_argument('--train_weight_path',type=str,default = None)

### Training Information File
parser.add_argument('--train_info_path',type=str,default = None)

### Genotype File
parser.add_argument('--genofile', type=str, default=None)
parser.add_argument('--genofile_header', type=str, default=None)

### Prediction SampleIDs
parser.add_argument('--sampleID',type=str,default=None)

### Specified Input Genotype File Type (vcf or dosages)
parser.add_argument('--geno', choices=['vcf', 'dosages'], type=str,default = None)

### GT or DS
parser.add_argument('--Format',choices=['GT', 'DS'],type=str,default=None)

### Window Size around Gene Boundary
parser.add_argument('--window',type=int,default=None)

### Threshold of Difference between Minor Allele Frequency (maf) Calculated in Training and Prediction
parser.add_argument('--maf_diff',type=float,default=None)

### Number of Thread
parser.add_argument('--thread',type=int,default = None)

### Output Dir
parser.add_argument('--out_prefix',type=str,default=None)

args = parser.parse_args()

############################################################################################
### Check Input Commands
print(args.train_weight_path)
print(args.train_info_path)
print("Using "+args.model+" in Training Part.")
print("Prediction sampleID:"+args.sampleID)
print("Chromosome Number:"+str(args.chr_num))
print("Prediction Genotype File:"+args.genofile)
print(args.genofile_header)

if args.geno=='vcf':
    print("Using "+args.Format+" Format for Prediction.")
elif args.geno=='dosages':
    print("Using DS Format for Prediction.")

print("Difference of MAF Threshold for Dropping Prediction snps:"+str(args.maf_diff))
print("window="+str(args.window))
print("Number of Thread:"+str(args.thread))
print("Output Dir:"+args.out_prefix)
############################################################################################
### Training Weight
Weight=pd.read_csv(args.train_weight_path,sep='\t')
Weight['CHROM']=Weight['CHROM'].astype('int')
Weight['POS']=Weight['POS'].astype('int')

### Create snpID (CHROM:POS:REF:ALT)
Weight['snpID']=(Weight['CHROM'].astype('str')+':'+Weight['POS'].astype('str')
                 +':'+Weight.REF+':'+Weight.ALT)

### Create Reversed snpID (CHROM:POS:ALT:REF)
Weight['snpID_Reversed'] = (Weight['CHROM'].astype('str')+':'+Weight['POS'].astype('str')
                            +':'+Weight.ALT+':'+Weight.REF)

if len(Weight.ES)==0:
    raise SystemExit('No Training Parameters. Please Check Input Files.')

### Training Information
Train_Info=pd.read_csv(args.train_info_path,sep='\t')

### Read in Header for Genotype File
genofile_header=pd.read_csv(args.genofile_header,sep='\t').rename(columns={'#CHROM':'CHROM'})
genofile_header=np.array(tuple(genofile_header))

### SampleIDs
sampleID=pd.read_csv(args.sampleID,sep='\t',header=None)
sampleID=np.array(sampleID).ravel()

### Gene with Weight
TargetID=unique(np.array(Weight.TargetID))

### Initializing Header of Output Files
pd.DataFrame(columns=np.hstack((['CHROM','GeneStart','GeneEnd','TargetID','GeneName'],
                                sampleID))).to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_'+args.model+'_GReX_prediction.txt',
                                                   sep='\t',index=None,header=True,mode='a')

#############################################################################################
### Thread Process
def thread_process(num):
    Info_temp=Train_Info >> mask(Train_Info.TargetID==TargetID[num])
    Info_temp=Info_temp.reset_index(drop=True)
    
    ### Gene Boundary
    start = max(int(Info_temp.GeneStart)-args.window,0)
    end = max(int(Info_temp.GeneEnd)+args.window,0)
    
    ### Select Corresponding Genotype File by Tabix
    chr_process=subprocess.Popen(["tabix"+" "+args.genofile+" "+str(args.chr_num)+":"+str(start)+"-"+str(end)],
                                  shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    
    out=chr_process.communicate()[0]
    
    if len(out)==0:
        print("No Testing Data for Gene:"+TargetID[num])
    else:
        print("Running Prediction for Gene:"+TargetID[num])
        ### Decode subprocess output with 'utf-8'
        Chr_temp=pd.read_csv(StringIO(out.decode('utf-8')),sep='\t',header=None,low_memory=False)
        Chr_temp.columns = np.array(tuple(genofile_header))
        Chr_temp=Chr_temp.reset_index(drop=True)
        
        ### Reform Genotype Data
        Chrom = CHR_Reform('PRED', Chr_temp, args.geno, sampleID, args.Format).main()
        
        ### Dealing with Training Weight
        beta_temp = (Weight >> mask(Weight.TargetID==TargetID[num]) 
                     >> select(Weight.snpID,Weight.snpID_Reversed,Weight.ES,Weight.MAF))
        beta_temp = beta_temp.drop_duplicates(['snpID'],keep='first')
        
        overlapID = np.intersect1d(np.array(beta_temp.snpID),np.array(Chr_temp.snpID))
        overlapID_Reversed = np.intersect1d(np.array(beta_temp.snpID_Reversed),np.array(Chr_temp.snpID))
        
        if (len(overlapID) == 0) & (len(overlapID_Reversed)==0):
            print("No Overlapping snpID Provided between Genotype and Training Weight File.")
        
        else:
            if len(overlapID)>=len(overlapID_Reversed):
                print("Number of Overlapping snpID:"+str(len(overlapID)))
                beta_temp = beta_temp >> drop(beta_temp['snpID_Reversed'])
            else:
                print("Number of Overlapping snpID:"+str(len(overlapID_Reversed)))
                beta_temp = (beta_temp >> drop(beta_temp['snpID'])).rename(columns={'snpID_Reversed':'snpID'})
                Chrom['MAF_test']=1-Chrom['MAF_test'].astype('float')
                
            Chrom = Chrom >> select(Chrom.snpID, Chrom.MAF_test, Chrom[sampleID])
            
            ### Store Prediction Result
            pred = pd.DataFrame()
            pred['sampleID'] = np.array(sampleID).ravel()
            
            Pred = (Chrom
                    >> mask(Chrom.snpID.isin(overlapID))).merge((beta_temp
                                                                 >> mask(beta_temp.snpID.isin(overlapID))),
                                                                left_on='snpID',
                                                                right_on='snpID',
                                                                how='outer')
            
            Pred['diff'] = abs(Pred['MAF'].astype('float')-Pred['MAF_test'].astype('float'))
            Pred = Pred >> mask(Pred['diff']<=args.maf_diff) >> drop(Pred[['MAF','MAF_test','diff']])
            print("Overall ID used:"+str(len(Pred.snpID)))
            
            if len(Pred.snpID)==0:
                print("Prediction Stop.")
                print("No snps Satisfied Condition Below:")
                print("Differences of maf between Training and Prediction less than "+str(args.maf_diff))
            else:
                testX_temp=Pred.T
                testX_temp.columns=testX_temp.loc['snpID']
                testX=pd.DataFrame(testX_temp.drop(['snpID','ES']),dtype='float')
                pred[TargetID[num]]=mat(testX)*(mat(Pred['ES'].astype('float')).T)
                
                pred=pred.T.reset_index(drop=True)
                pred.columns=pred.loc[0]
                pred=pred.drop([0])
                pred['TargetID']=TargetID[num]
                
                out=(Info_temp[['CHROM','GeneStart','GeneEnd',
                                'TargetID','GeneName']].merge(pred,left_on='TargetID',
                                                              right_on='TargetID',how='outer'))

                out.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_'+args.model+'_GReX_prediction.txt',
                           sep='\t',index=None,header=None,mode='a')

########################################################################################################
# Start Thread
if (args.thread < int(len(TargetID)/100) | args.thread > len(TargetID)):
    args.thread = (int(len(EXP)/100)+1)*100

pool = multiprocessing.Pool(args.thread)

pool.map(thread_process,[num for num in range(len(TargetID))])

pool.close()

pool.join()

########################################################################################################
### Time Calculation
end_time=time()
print("Running Time:"+str(round((end_time-start_time)/60,2))+" minutes.")









