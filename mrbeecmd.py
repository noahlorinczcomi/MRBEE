print('Getting things ready')
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas
import numpy
from statsmodels import regression
from numpy import matrix
from scipy import stats
import gzip
import argparse
from mrbeecmdfunctions import *

parser=argparse.ArgumentParser(prog='MRBEE using GWAS summary statistics',
                               description='This program performs univariable and multivariable MR with MRBEE using GWAS summary statistics',
                              allow_abbrev=False)

# required
parser.add_argument('-exposureGWAS',action='store',type=str,help='(Required) Comma-separated filepaths (w/ .gz extensions) to all exposure GWAS data sets. You may simply put one filepath for univariable MR')
parser.add_argument('-outcomeGWAS',action='store',type=str,help='(Required) Filepath (w/ .gz extension) to phenotype GWAS data set')
parser.add_argument('-IVs',action='store',type=str,help='(Required) Path to file (w/ extension) containing unique SNP identifiers for IVs you wish to use in MR. These unique SNP identifiers must be in the same format as those for each exposure and the outcome.')

parser.add_argument('-exposureSNP',action='store',type=str,default='',help='Comma-separated list of unique SNP identifier column names in each of the exposure GWAS data sets. Order must match order of GWAS data set filepaths.')
parser.add_argument('-outcomeSNP',action='store',type=str,default='',help='Column name of unique SNP identifier in the outcome GWAS data set.')

parser.add_argument('-exposureBETA',action='store',type=str,default='',help='Comma-separated list of BETA (association estimate/effect size) column names in each of the exposure GWAS data sets. Order must match order of GWAS data set filepaths. Assumed to already be in standardized scale unless the --stdZ True flag is included.')
parser.add_argument('-outcomeBETA',action='store',type=str,default='BETA',help='BETA (association estimate/effect size) column name in the outcome GWAS data set. Assumed to already be in standardized scale unless the --stdZ True flag is included.')

parser.add_argument('-exposureSE',action='store',type=str,default='',help='Comma-separated list of standard error column names in each of the exposure GWAS data sets. Order must match order of GWAS data set filepaths.  Assumed to already be in standardized scale unless the --stdZ True flag is included.')
parser.add_argument('-outcomeSE',action='store',type=str,default='SE',help='Standard error column name in the outcome GWAS data sets. Assumed to already be in standardized scale unless the --stdZ True flag is included.')

parser.add_argument('-exposureEffectAllele',action='store',type=str,default='',help='Comma-separated list of effect allele column names in each of the exposure GWAS data sets. Order must match order of GWAS data set filepaths.')
parser.add_argument('-outcomeEffectAllele',action='store',type=str,default='ALT',help='(Required) Standard error column name in the outcome GWAS data sets. Assumed to already be in standardized scale unless the --stdZ True flag is included.')

parser.add_argument('-stdZ',action='store',type=str,default='True',help='Should standardization on BETA and SE be performed such that new standardized effect sizes are equal to BETA/SE and their standard errors are 1? If yes, put True.')

parser.add_argument('-exposureNames',action='store',type=str,default='',help='Comma-separated list of exposure names, ordered according to the order of specified exposure filepaths. If you do not specify this, the program will assign names such as exposure 1, exposure 2, etc.')

parser.add_argument('-genomewideSpleio',action='store',type=str,default='False',help='Do you want to perform genome-wide horizontal pleiotropy testing? Put "True" if yes, otherwise put "False", which is  the default')

parser.add_argument('-out',action='store',type=str,default='results',help='Filepath (w/o extension) to which results should be written. Default is "results"')


args=parser.parse_args()

## checking on defaults
p=len(args.exposureGWAS.split(','))

args.exposureSNP=usedefault(args.exposureSNP,('SNP,'*p)[:-1]) # exposure SNPs
args.exposureBETA=usedefault(args.exposureBETA,('BETA,'*p)[:-1]) # exposure BETAs
args.exposureSE=usedefault(args.exposureSE,('SE,'*p)[:-1]) # exposure SEs

# this is the simple MRBEE command line tool
# steps
# 1) load IVs
# 2) load all GWAS data and subset to only IVs
# 3) merge subsetted GWAS data
# 4) harmonise alleles
# 5) calculate bias-correction terms
# 6) subset working data to only IVs
# 7) perform MR
# 8) save results

# 1) load IVs
ivs=pandas.read_csv(args.IVs, engine='python', header=0).values.squeeze()

# 2), 3) load all GWAS data (do not subset yet)
uc=[args.outcomeSNP,args.outcomeBETA,args.outcomeSE,args.outcomeEffectAllele] # outcome column names
# checking if all column names user gave are actually in data ...
# code here ...
print('Reading in outcome GWAS data')
outcomedata=pandas.read_csv(args.outcomeGWAS,compression='gzip',header=0,delim_whitespace=True,usecols=uc)
outcomedata=outcomedata.rename(columns={uc[0]: 'SNP', uc[1]: 'BETAy', uc[2]: 'SEy', uc[3]: 'EAy'}) # rename to names I'll know later

efps=args.exposureGWAS.split(','); p=len(efps) # list of ordered exposure filepaths
snps=args.exposureSNP.split(',') # list of exposure SNP column names
betas=args.exposureBETA.split(',') # ...
ses=args.exposureSE.split(',') # ...
eas=args.exposureEffectAllele.split(',') # ...
bdf=outcomedata.copy()
for _ in range(0,len(efps)):
    print('Reading in GWAS data for exposure '+str(_+1))
    # checking if all column names user gave are actually in data ...
    # code here ...
    data=pandas.read_csv(efps[_],compression='gzip',header=0,delim_whitespace=True,usecols=[snps[_],betas[_],ses[_],eas[_]])
    data=data.rename(columns={snps[_]: 'SNP',betas[_]: 'BETAx'+str(_), ses[_]: 'SEx'+str(_), eas[_]: 'EAx'+str(_)}) # rename to names I'll know later
    bdf=pandas.merge(bdf,data,how='inner',left_on='SNP',right_on='SNP')
    if _==0:
        print('1 exposure down, '+str(p-1)+' to go')
    else:
        print(str(_+1)+' exposures down, '+str(p-_-1)+' to go')

# 4) harmonise all alleles (effect sizes)
print('Harmonizing alleles')
for _ in range(0,p):
    betacn='BETAx'+str(_)
    eacn='EAx'+str(_)
    mask=bdf[eacn].str.upper()!=bdf['EAy'].str.upper() # find where outcome effect allele does not match this exposure's effect allele
    bdf.loc[mask,betacn]=(-1*bdf.loc[mask,betacn]) # harmonise to outcome effect allele

# 5) bias-correction correlation matrix
print('Calculating bias-correction terms')
betanames=['BETAy']
senames=['SEy']
for _ in range(0,p):
    betanames.append('BETAx'+str(_))
    senames.append('SEx'+str(_))

# Note that I changed the code below slightly - I don't want users subsampling bc I noticed the results can vary
if bdf.shape[0]>float('inf'): # if data is very large, subset to random sample of 500k
    ss=bdf.sample(n=int(5e5),axis=0)
    inds=numpy.zeros((int(5e5),p+1),dtype='bool')
    for _ in range(0,p+1):
        chi=(ss[betanames[_]]/ss[senames[_]])**2
        mask=chi<3.841459 # these are rows to use
        inds[mask,_]=True
    dokeep=numpy.sum(inds,axis=1)>0
    cc=ss[dokeep][betanames].corr()
else: # else, use all data
    inds=numpy.zeros((bdf.shape[0],p+1),dtype='bool')
    for _ in range(0,p+1):
        chi=(bdf[betanames[_]]/bdf[senames[_]])**2
        mask=chi<3.841459 # these are rows to use
        inds[mask,_]=True
    dokeep=numpy.sum(inds,axis=1)>0
    cc=bdf[dokeep][betanames].corr()

# 6) now can subset to only IVs
if args.genomewideSpleio.title()=='True': 
    bdfog=bdf.copy() # save a copy only if user also wants to perform genome-wide Spleio testing (if they don't not saving a copy saves space)sssssssssss
bdf=bdf[bdf['SNP'].isin(ivs)]

# 7) perform MR
print('Performing MR')
m=bdf.shape[0]
ebeta=['BETAx'+str(_) for _ in range(0,p)]
ese=['SEx'+str(_) for _ in range(0,p)]
bx=bdf[ebeta].values.reshape((m,p))
bxse=bdf[ese].values.reshape((m,p))
by=bdf['BETAy'].values.reshape((m,1))
byse=bdf['SEy'].values.reshape((m,1))

if args.stdZ.title()=='True': # if user wants the Z-score standardization
    bx=bx/bxse
    bxse=bxse/bxse
    by=by/byse
    byse=byse/byse

# bias-correction Sigma
meanses=numpy.median(numpy.column_stack((byse,bxse)),axis=0)
D=numpy.diag(meanses.squeeze())
Sigma=D@cc@D # cc is a pandas object so can use .iloc
VV=numpy.array(Sigma.iloc[0,0]).reshape((1,1))
UU=Sigma.iloc[1:,1:].values
UV=Sigma.iloc[1:,0].values.reshape((p,1))

I=numpy.eye(m)
print('Using '+str(m)+' candidate SNPs in MR')
est,V,outliers,ign,kiter=imrbee(bx,by,UU,UV,VV,I,I,PleioPThreshold=0.05/m**0.5)
print(str(len(outliers))+' SNPs removed from MR because of horizontal pleiotropy, '+str(m-len(outliers))+' SNPs remain')

cn=['intercept']
if len(args.exposureNames)>0:
    cn.append(args.exposureNames.split(','))
elif (len(args.exposureNames)!=p) | (len(args.exposureNames)==0):
    myown=['Exposure'+str(_) for _ in range(0,p)]
    cn.append(myown)

cn=flatten_list(cn)
est=est.squeeze()
ses=(numpy.diag(V))**0.5
P=[1-stats.chi2.cdf((est[_]/ses[_])**2,1) for _ in range(0,p+1)]
df=pandas.DataFrame({'Exposure': cn, 'Est': est, 'SE': ses, 'P': P})
print(df)
mask=numpy.ones((m,),dtype='bool'); mask[outliers]=False # mask for if IV was used by MRBEE and not excluded by the pleiotropy test
UU_=numpy.column_stack(([0]*p,UU)); UU_=numpy.row_stack(([0]*(p+1),UU_))
UV_=numpy.row_stack((0,UV))
finalsp=SpleioP(numpy.column_stack(([1]*m,bx)),by,UU_,UV_,VV,est,V) # final Spleio P-values
# does the user also want to perform genome-wide horizontal pleiotropy testing?
if args.genomewideSpleio.title()=='True':
    print('Performing genome-wide horizontal pleiotropy testing using Spleio')
    BX=bdfog[ebeta].values
    BXSE=bdfog[ese].values
    BY=bdfog['BETAy'].values
    BYSE=bdfog['SEy'].values
    if args.stdZ.title()=='True':
        BX=BX/BXSE
        BY=BY/BYSE
        BXSE=numpy.ones((M,p))
        BYSE=numpy.ones((M,1))
    RESDF=genomePleio(BX,BY,BXSE,BYSE,cc,est,V)

newcn1=['by']
newcn2=['byse']
for _ in range(0,p):
    newcn1.append('bx'+str(_+1))
    newcn2.append('bxse'+str(_+1))

newcn=flatten_list(['SNP',newcn1,newcn2,'notDroppedBySpleio','finalSpleioPvalues'])
meta=pandas.DataFrame(numpy.column_stack((bdf['SNP'].values,by,bx,byse,bxse,mask,finalsp)),columns=newcn)

# save causal estimates (`df`), MR data (`meta`), genome-wide pleiotropy testing (`RESDF`)
print('Saving results to {}_<causal_estimates,MR_data,genomewide_Spleio>.csv'.format(args.out))
df.to_csv(args.out+'_causal_estimates.csv')
meta.to_csv(args.out+'_MR_data.csv')
if args.genomewideSpleio.title()=='True':
    RESDF.to_csv(args.out+'_genomewide_Spleio.csv')


