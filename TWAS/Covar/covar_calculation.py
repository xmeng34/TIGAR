#!/usr/bin/env python

##############################################################################################
# Import packages needed
import argparse
from time import *
import subprocess
from subprocess import *
import pandas as pd
import numpy as np
from numpy import *
from dfply import *
import io
from io import StringIO
from io import *

import multiprocessing
###############################################################################################
# Hide Warnings from Python
warnings.filterwarnings("ignore")
###############################################################################################
### Time Calculation
start_time=time()
###############################################################################################
# Reform vcf file
### input each sample genotype
### For GT Format:
###  code '0|0' or '0/0' as 0
###  code ('0|1' or '1|0')  or ('0/1' or '1/0') as 1
###  code '1|1' or '1/1' as 2
###  code '.|.' or './.' as nan(missing)

### For DS Format:
### code '.' as nan(missing)
def geno_reform(data,Format):
    if Format=='GT':
        data[(data=='0|0')|(data=='0/0')]=0
        data[(data=='1|0')|(data=='1/0')|(data=='0|1')|(data=='0/1')]=1
        data[(data=='1|1')|(data=='1/1')]=2
        data[(data=='.|.')|(data=='./.')]=nan
    elif Format=='DS':
        data[(data==".")]=nan
    return data

##########################################################################################
### For vcf input
### Split input dataframe by Format. ex, '0|0:0.128'
### Input:
### 1. data:The first nine columns fixed
### 2. Format: GT or DS

### Output:
###  The First six columns of output dataframe should be:
###    1) CHROM
###    2) POS
###    3) ID (i.e. rsID)
###    4) REF
###    5) ALT
###    6) snpID (CHROM:POS:REF:ALT)
###    7) p_HWE:p-value for Hardy Weinberg Equilibrium exact test
###    8) MAF: Minor Allele Frequency (range from 0~1)
###    9) Samples gene variance splited by Format (GT or DS)

def CHR_Reform(data, geno, Format, maf):
    if geno=='vcf':
        sampleID = data.columns[9:]
    elif geno=='dosages':
        sampleID = data.columns[5:]
    
    data['snpID']=(data['CHROM'].astype('str')+":"+data['POS'].astype('str')
                   +":"+data.REF+":"+data.ALT)
        
    CHR = data >> select(data[['CHROM','POS','ID','REF','ALT','snpID']],data[sampleID])

    CHR=CHR.drop_duplicates(['snpID'],keep='first')
    
    if geno=='vcf':
        indicate=data.FORMAT[0].split(":").index(Format)
        CHR[sampleID]=CHR[sampleID].applymap(lambda x:x.split(":")[indicate])
    
    CHR[sampleID]=CHR[sampleID].apply(lambda x:geno_reform(x,Format),axis=0)

    ### Calculate MAF by SNPs(range from 0-1)
    temp=pd.DataFrame((CHR >> select(CHR[sampleID])),dtype=np.float)
    CHR['MAF']=temp.apply(lambda x:sum(x)/(2*len(x.dropna())),axis=1)

    ### Dealing with NaN
    CHR[np.hstack(([sampleID,'MAF']))] = CHR[np.hstack(([sampleID,'MAF']))].apply(lambda x:x.fillna(2*x.MAF),axis=1)
    
    return (CHR>>mask(CHR.MAF>maf)),sampleID

######################################################################################
### variable needed
parser = argparse.ArgumentParser(description='manual to this script')

### Chromosome Block Information
parser.add_argument('--block',type=str,default = None)

### Genotype File
parser.add_argument('--genofile',type=str,default = None)
parser.add_argument('--genofile_header',type=str, default=None)

### Specified Input Genotype File Type (vcf or dosages)
parser.add_argument('--geno',type=str, choices=['vcf', 'dosages'], default = None)

### Chromosome Number
parser.add_argument('--chr_num',type=int,default = None)

### 'DS' or 'GT'
parser.add_argument('--Format',choices=['GT', 'DS'],type=str,default = None)

### Minor Allele Frequency (maf, Range from 0-1) Threshold 
### for Seleting Genotype Data to Calculate Covariance Matrix
parser.add_argument('--maf',type=float,default = None)

### Number of Thread
parser.add_argument('--thread',type=int,default = None)

### Output Dir
parser.add_argument('--out_prefix',type=str,default = None)

args = parser.parse_args()

######################################################################################
### variable checking
print("Block information:"+args.block)
print("Genotype File:"+args.genofile)
print(args.genofile_header)

if args.geno=='vcf':
    print("Using "+args.Format+" Format for association study.")
elif args.geno=='dosages':
    print("Using DS format for association study.")

print("Chromosome number:"+str(args.chr_num))
print("Number of thread:"+str(args.thread))
print("Threshold of maf value:"+str(args.maf))
print("Output dir:"+args.out_prefix)

#######################################################################################
### Read in Block Information
Block = pd.read_csv(args.block,sep='\t')

### Read in Header for Genotype File
genofile_header=pd.read_csv(args.genofile_header,sep='\t').rename(columns={'#CHROM':'CHROM'})

### Initializing Header of Output Files
pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT','COV']).to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_reference_cov.txt',
                                                                        sep='\t',index=None,header=True)
#########################################################################################
def thread_process(num):
    block_temp = Block.loc[num]
    
    ### Select Corresponding Genotype File by Tabix
    chr_process=subprocess.Popen(["tabix"+" "+args.genofile+" "+str(args.chr_num)+":"+str(block_temp.Start)+"-"+str(block_temp.End)],
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
    out=chr_process.communicate()[0]

    if len(out)==0:
        print("No corresponding genotype data in block:"+str(block_temp.Start)+'--'+str(block_temp.End))
    else:
        ### Decode subprocess output with 'utf-8'
        Chr_temp = pd.read_csv(StringIO(out.decode('utf-8')),sep='\t',low_memory=False)
        Chr_temp.columns = np.array(tuple(genofile_header))
        Chr_temp = Chr_temp.reset_index(drop=True)

        Chrom, sampleID = CHR_Reform(Chr_temp, args.geno, args.Format, args.maf)
        Chrom = Chrom.sort_values(by='POS').reset_index(drop=True)

        ### Calculating Correlation Coefficient Matrix
        mcovar = cov(pd.DataFrame(Chrom[sampleID],dtype='float'))
        
        for i in range(len(Chrom)):
            covar_info=np.append([Chrom.loc[i][0:5].ravel()],','.join(mcovar[i,i:].astype('str')))
            pd.DataFrame(covar_info).T.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_reference_cov.txt',
                                              sep='\t',index=None,header=None,mode='a')

########################################################################################################
### Start Thread
pool = multiprocessing.Pool(args.thread)

pool.map(thread_process,[num for num in range(len(Block))])

pool.close()
pool.join()

#########################################################################################################
### Time Calculation
end_time=time()
print("Running Time:"+str(round((end_time-start_time)/60,2))+" minutes.")















