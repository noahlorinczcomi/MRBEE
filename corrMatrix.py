
import argparse
import pandas
from scipy.stats import norm

parser = argparse.ArgumentParser(
                    prog='Calculation of Correlation',
                    description='This program calculates the correlation between GWAS estimates for a pair of phenotypes')

# now trying to allow for any number of data sets and column names to be passed
parser.add_argument('-data',action='store',type=str,help='(Required) comma-separated list of filepaths to all GWAS sumstats')
parser.add_argument('-snp',action='store',type=str,help='(Required) comma-separated list of SNP column names for all GWAS sumstats')
parser.add_argument('-beta',action='store',type=str,help='(Required) comma-separated list of names of BETA columns in GWAS sumstats')
parser.add_argument('-se',action='store',type=str,help='(Required) comma-separated list of names of SE columns in GWAS sumstats')
parser.add_argument('-pt',action='store',default=0.05,type=float,help='(Required) P-value threshold. Only SNPs with all GWAS estimates P>this threshold will be considered')
parser.add_argument('-names',action='store',type=str,help='(Required) Comma-separated list of row, column names to assign to the correlation matrix')
parser.add_argument('-out',action='store',type=str,help='(Required) Directory (file location w/o extension) of location to write out correlation matrix')

args=parser.parse_args()

print("NOTE: -snp, -beta, -se, and -names flag declarations must be in the corresponding order of -data declarations")

fps=args.data; dsplitted=fps.split(",")
snps=args.snp; snpsplitted=snps.split(",")
bs=args.beta; bsplitted=bs.split(",")
ss=args.se; ssplitted=ss.split(",")

datadict={};names=[]
for i in range(0,len(dsplitted)):
    nn="data{}".format(i); names.append(nn)
    uc=[snpsplitted[i], bsplitted[i], ssplitted[i]] # don't need all columns (will speed later things up)
    dat=pandas.read_table(dsplitted[i],sep=None,engine='python',usecols=uc) # will auto find the delimiter
    dat["z{}".format(i)]=dat[bsplitted[i]]/dat[ssplitted[i]]
    # join to earlier dats
    if i>0:
        prevkey="data{}".format(i-1)
        prevdat=datadict[prevkey]
        prevsnp=snpsplitted[i-1]; currsnp=snpsplitted[i]
        merged=pandas.merge(left=prevdat, right=dat, left_on=prevsnp, right_on=currsnp)
        # need also to assign current dict entry as merged data so it can be merged again next
        datadict[nn]=merged
        names.append(nn)
        del dat
    else:
        datadict[nn]=dat
        names.append(nn)
        del dat

# now calculate correlation
allmerged=datadict[names[len(names)-1]]
# make a matrix of just Z-statistics, preserving original order
zn=["z{}".format(i) for i in range(0, len(dsplitted))]
Z=allmerged[zn]
# removing rows where any |A|>pt
Zcut=norm.ppf(1-args.pt/2)
def dropRows(Zmat,Zcut):
    isKeep=[]
    for i in range(0,Zmat.shape[0]):
        if all((abs(Zmat.loc[i,:])>Zcut)==False)==True:
            isKeep.append(i)
    return isKeep

isKeep=dropRows(Z,Zcut)
print("the program is running")
R=Z.loc[isKeep,:].corr();
# changing row, column names of correlation matrix
nps=args.names; nsplitted=nps.split(",")
R.index=nsplitted; R=R.set_axis(list(nsplitted), axis=1, inplace=False)
R.to_csv(args.out, index=True,header=True,sep=' ')
print(R)
print('{} SNPs used in correlation matrix estimation'.format(Z.loc[isKeep,:].shape[0]))

# testing
"""
fps="C:/Users/njl96/Downloads/d1b.csv,C:/Users/njl96/Downloads/d2.csv,C:/Users/njl96/Downloads/d3.csv"
fps="/Users/Noah/Downloads/d1.csv.gz,/Users/Noah/Downloads/d2.csv.gz,/Users/Noah/Downloads/d3.csv.gz"
fps="/Users/Noah/Downloads/d1b.txt,/Users/Noah/Downloads/d2b.txt,/Users/Noah/Downloads/d3b.txt"
dsplitted=fps.split(",")
snpsplitted="X,X,X".split(",")
bsplitted="b1,b2,b3".split(",")
ssplitted="s1,s2,s3".split(",")
"""





