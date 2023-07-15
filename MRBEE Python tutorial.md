# MRBEE using Python
Please note that this tutorial has been successfully run using Python v3.7.0 and requires the following modules: `pandas`, `numpy`, `scipy`, `statsmodels`, `gzip`. To retrieve the data used in this tutorial, execute the following command in a linux system:
```unix
wget -O PyTutorialData.tar.gz "https://www.dropbox.com/s/2fbnonxt3xds3to/PyTutorialData.tar.gz?dl=0"
tar -xf PyTutorialData.tar.gz
wget https://raw.githubusercontent.com/noahlorinczcomi/MRBEE/main/mrbeecmdfunctions.py
wget https://raw.githubusercontent.com/noahlorinczcomi/MRBEE/main/mrbeecmd.py
```
The `mrbeecmd.py` file is the MRBEE command line tool and is not necessary to complete this tutorial. Please see Lorincz-Comi & Yang et al. (2023) for a demonstration of how to use the `mrbeecmd.py` command line tool.
This tutorial demonstrates how to use the MRBEE Python module and is divided into four sections:
1. Loading and harmonising GWAS data
2. Calculating bias-correction terms
3. Performing MR using MRBEE
4. Performing genome-wide horizontal pleiotropy testing
Before we complete any of these steps, we must first source the MRBEE functions from the `mrbeecmdfunctions.py` file, which we assume to be located in the current working directory.
```python
from mrbeecmdfunctions import *
```
# Loading and harmonising GWAS data
In this section, we assume that all GWAS summary statistic data sets are loaded in the current working directory and each have the '.gz' file extension. The gz file compression can simply be added to a non-compressed file using 'gzip file' on a Linux system.
We have three data sets: `cad_gwas.txt.gz`, `hba1c_gwas.txt.gz`, `uricacid_gwas.txt.gz`, where the first corresponds to coronary artery disease in Europeans from the CARDIoGRAMplusC4D cohort and the latter two respectively to HbA1c and uric acid in Europeans in the UK Biobank. We first load the CAD data, then we load the list of IVs that we eventually want to try using in MR. These IVs are in the `IVs.txt` file in the same directory, and after loading them we will convert them to a vector. This vector will be a numpy array of rsIDs.
```python
uc=['markername','beta','se_dgc','effect_allele'] # defining names of only select coolumns I want to load to save memory
outcomedata=pandas.read_csv('cad_gwas.txt.gz',compression='gzip',header=0,delim_whitespace=True,usecols=uc)
outcomedata=outcomedata.rename(columns={uc[1]: 'BETAy', uc[2]: 'SEy', uc[3]: 'EAy'}) # renaming to more informative names
IVs=pandas.read_csv('IVs.txt', engine='python', header=0).values.squeeze() # len(IVs)==1520, type(IVs)==<class 'numpy.ndarray'>
```
Now that the outcome GWAS data is loaded in the `outcomedata` object, we can load the two exposure GWAS data sets and sequentially merge them with `outcomedata` by the `markername` column in `outcomedata`. The `markername` column in `outcomedata` is the unique SNP identifier column of rsIDs. The two exposure GWAS data sets also have columns containing rsIDs, though the names of these columns are `rsid` for both exposures and not `markername`. As a result, when each exposure GWAS is loaded, we will change the `rsid` column name to `markername`.
```python
ucx=['rsid','beta','se','alt'] # both exposure GWAS data sets have the same column names
datax1=pandas.read_csv('hba1c_gwas.txt.gz',compression='gzip',header=0,delim_whitespace=True,usecols=ucx)
datax1=datax1.rename(columns={'rsid': 'markername', 'beta': 'BETAx1', 'se': 'SEx1', 'alt': 'EAx1'}) # provide informative names
bdf=pandas.merge(outcomedata,datax1,how='inner',left_on='markername',right_on='markername') # merge HbA1c and CAD GWAS
datax2=pandas.read_csv('uricacid_gwas.txt.gz',compression='gzip',header=0,delim_whitespace=True,usecols=ucx)
datax2=datax2.rename(columns={'rsid': 'markername', 'beta': 'BETAx2', 'se': 'SEx2', 'alt': 'EAx2'}) # provide informative names
bdf=pandas.merge(bdf,datax2,how='inner',left_on='markername',right_on='markername') # merge all exposure and outcome GWAS
del outcomedata, datax1, datax2 # no longer need these
```
The object `bdf` contains all GWAS summary statistics for HbA1c, uric acid, and CAD merged together by rsID (i.e., the `markername` column in `bdf`). The next step is to harmonise the exposure effect alleles to the outcome effect alleles. This process ensures that effect size estimates for the same SNP across phenotypes each correspond to increases in counts of the same allele.
```python
mask1=bdf['EAx1'].str.upper()!=bdf['EAy'].str.upper() # True where HbA1c and CAD effect alleles match, False elsewhwere
mask2=bdf['EAx2'].str.upper()!=bdf['EAy'].str.upper() # True where uric acid and CAD effect alleles match, False elsewhwere
bdf.loc[mask1,'BETAx1']=(-1*bdf.loc[mask1,'BETAx1']) # if not harmonised, change sign of HbA1c effect size
bdf.loc[mask2,'BETAx2']=(-1*bdf.loc[mask2,'BETAx2']) # if not harmonised, change sign of uric acid effect size
```
# Calculating bias-correction terms
Now we can calculate the correlation matrix **R** of measurement errors for each GWAS. This is the primary component of the bias-correction terms. Eventually, we will transform **R** to be a variance-covariance matrix using the standard errors of the exposure and outcome GWAS effect size estimates. These are in the `SEx1`, `SEx2`, and `SEy` columns for HbA1c, uric acid, and CAD, respectively.
Calculation of **R** is done using the Pearson correlation estimate and only SNPs with P>0.05for HbA1c, uric acid, and CAD in GWAS (Zhu et al., 2015; AJHG). The first step in this process is finding those SNPs with P>0.05. An example of how we can do this is below.
```python
inds=numpy.zeros((bdf.shape[0],3),dtype='bool')
betanames=['BETAy','BETAx1','BETAx2'] # corresponds to CAD, HbA1c, uric acid
senames=['SEy','SEx1','SEx2']
# I am using a for loop to save space, looping over each phenotype
for _ in range(0,3):
    chi=(bdf[betanames[_]]/bdf[senames[_]])**2
    mask=chi<3.841459 # 3.841459 is the 95th quantile of the central chi-square distribution with 1 degree of freedom
    inds[mask,_]=True # <2 seconds for ~8M SNPs and 3 phenotypes

snpsToUse=numpy.sum(inds,axis=1)>0 # where True, these SNPs have P>0.05 for CAD, HbA1c, and uric acid in each of their original GWAS
```
Now we have the `snpsToUse` object which is a vector of True or False values for each SNP. If `snpsToUse` is True for a particular SNP, that means we will use the SNP to estimate **R**, otherwise we will not. Upon checking, we see that over 8.1 million SNPs have True in `snpsToUse`. The next step is to calculate the matrix of Pearson correlations between each phenotype, only using SNPs where `snpsToUse` is True, which is done like this:
```python
R=bdf[snpsToUse][['BETAy','BETAx1','BETAx2']].corr()
```
As mentioned, now that we have the correlation matrix **R**, we need to use the standard errors in the corresponding columns `SEy`, `SEx1`, and `SEx2` to convert **R** to **Σ**, the corresponding variance-covariance matrix. This means that we need to extract standard errors for CAD, HbA1c, and uric acid. In practice, we recommend doing this step for only the instrumental variables that you intend to use in MR. As such, in the next step we convert **R** to **Σ** only for the 1,520 IVs we want to try using in MRBEE. First, we subset `bdf` to only the set of IVs and we name the new subsetted object as `mrdf`.
```python
mrdf=bdf[bdf['markername'].isin(IVs)]
```
Now, we can extract effect sizes and standard errors from `mrdf` to calculate **Σ** and get ready to perform MR. We perform these extractions like this:
```python
m=mrdf.shape[0] # number of IVs
p=2 # number of exposures
bx=mrdf[['BETAx1','BETAx2']].values.reshape((m,p)) # mxp matrix of unstandardized exposure effect sizes
bxse=mrdf[['SEx1','SEx2']].values.reshape((m,p))   # mxp matrix of unstandardized exposure standard errors
by=mrdf['BETAy'].values.reshape((m,1))             # mx1 vector of unstandardized outcome effect sizes
byse=mrdf['SEy'].values.reshape((m,1))             # mx1 vector of unstandardized outcome standard errors
```
Following the methods in Lorincz-Comi et al. (2023), we will use a Z-score standardization on `bx`, `bxse`, `by`, and `byse`. This is done with the following commands:
```python
bx=bx/bxse     # these are now exposure Z-scores/test-statistics
bxse=bxse/bxse # to make all 1s
by=by/byse     # these are now outcome Z-scores/test-statistics
byse=byse/byse # to make all 1s
```
The final step is to actually calculate the variance-covariance matrix of measurement errors **Σ**. Recall that the correlation matrix of measurement errors **R** has values corresponding to CAD in the first row/column, HbA1c in the second row/column, and uric acid (UA) in the third row/column. The transformation from **R** to **Σ** we are about to perform will have to consider that. We will perform this transformation using the matrix **D**, a diagonal matrix of averaged standard errors across the IVs for each phenotype. *Note*: Since we used the Z-score standardization, this step is unnecessary since **R**=**Σ**. However, for illustrative purposes, we demonstrate how to perform this step with your own data.
```python
bars=numpy.sum(numpy.column_stack((byse,bxse)),axis=0)/m
D=numpy.diag(bars.squeeze())
Sigma=D@R@D
```
We will now split `Sigma` (representing **Σ**) into three parts: (i) `SigmaUU` corresponding to measurement error covariance between the exposures, (ii) `SigmaUV` corresponding to measurement error covariance between the exposures and CAD, (iii) `SigmaVV` corresponding to measurement error variance of the outcome. Recall the row/column positions of Sigma/**Σ** occupied by CAD (1st,1st), HbA1c (2nd,2nd), and uric acid (3rd,3rd).
```python
SigmaVV=numpy.array(Sigma.iloc[0,0]).reshape((1,1))
SigmaUU=Sigma.iloc[1:,1:].values
SigmaUV=Sigma.iloc[1:,0].values.reshape((p,1))
```
# Performing MR using MRBEE
Now we have everything we need to perform MR using MRBEE. Currently, this software supports the MRBEE-IMRP method (see Lorincz-Comi, Yang, & Zhu, 2023), but future development will include the MRBEE-Median, MRBEE-MLqe, MRBEE-Mix, and MRBEE-IPOD methods. Each are almost identical in causal estimation but differ only in how horizontal pleiotropy is addressed. As such, the `imrbee` function loaded from the `mrbeecmdfunctions.py` file performs MRBEE-IMRP. MRBEE-IMRP requires specifying a parameter representing the P-value threshold below which specific IVs with evidence of horizontal pleiotropy using the $S_\text{pleio}$ statistic will be excluded from causal estimation.
Since we are using an independent set of IVs that are not in linkage disequilibrium with each other, we can use a Bonferroni-corrected P-value threshold, which in this case is $0.05/1445\approx 3.5\times 10^{-5}$ In practice, it may be more desirable to use a lower threshold, so in this case we will use $0.05/\sqrt{1445}\approx 1.3\times 10^{-3}$. This threshold will remove more horizontal pleiotropy than the $3.5\times 10^{-5}$ threshold. In contrast to the MRBEE functionality in the MRBEE R package, the `imrbee` function that we will use below to perform MRBEE-IMRP can accept an LD matrix of correlations between IVs. In our case, IVs are independent and so we will set this matrix to be equal to the identity matrix.
```python
I=numpy.eye(m)
est,V,outliers,ign,kiter=imrbee(bx,by,SigmaUU,SigmaUV,SigmaVV,I,I,PleioPThreshold=0.05/m**0.5)
```
The object `est` contains causal estimates for the 2 exposures and one leading intercept term. Note that it is generally advisable to include an intercept term in MVMR (see Lorincz-Comi et al., 2023) for consistent causal estimation. The object `V` is the $(2+1)\times (2+1)$ variance-covariance matrix corresponding to the causal estimates; `outliers` is a vector of index positions (from 0 to 1444, the number of IVs minus 1 because of Python's indexing scheme) of SNPs that had evidence of horizontal pleiotropy at the P-value threshold we selected; `ign` is the inverse of the LD matrix only for SNPs that were not outliers; `kiter` is the number of iterations the program required before converging, which in this case was 2. 
We can create a cleaner version of these results like this:
```python
est=est.squeeze()
ses=(numpy.diag(V))**0.5
P=[1-stats.chi2.cdf((est[_]/ses[_])**2,1) for _ in range(0,p+1)]
df=pandas.DataFrame({'Exposure': ['Intercept','HbA1c','Uric acid'], 'Est': est, 'SE': ses, 'P': P})
```
# Performing genome-wide horizontal pleiotropy testing
Genome-wide horizontal pleiotropy testing is a natural next step after causal estimation. Genome-wide horizontal pleiotropy testing can identify novel loci, pleiotropic loci, and loci whose association with CAD are completely mediated by HbA1c and/or uric acid (Lorincz-Comi et al., 2023; Zhu et al., 2022). Since this step requires some calculating of $\boldsymbol\Sigma_j$ for all of the *j*th SNPs genome-wide, it is better to use the memory-saving `genomePleio()` function that was initially sourced from the `mrbeecmdfunctions.py` file. First, we extract the genome-wide correlates of `bx`, `bxse`, `by`, and `byse` that we used earlier in MR. Now, we want to extract the same information but now genome-wide. **Importantly**, the same standardization must be applied in this step as was applied to the IVs used in MR. This is because $S_\text{pleio}$ uses the causal estimates and their corresponding variance-covariance matrix, which must be on the same scale as the data.
```python
M=bdf.shape[0] # number of SNPs genome-wide
BX=bdf[['BETAx1','BETAx2']].values
BXSE=bdf[['SEx1','SEx2']].values
BY=bdf['BETAy'].values.reshape((M,1))
BYSE=bdf['SEy'].values.reshape((M,1))
# standardization below
BX=BX/BXSE
BY=BY/BYSE
BXSE=numpy.ones(BXSE.shape)
BYSE=numpy.ones(BYSE.shape)
gwp=genomePleio(BX,BY,BXSE,BYSE,R,est,V) # R: correlation matrix of measurement errors; est: causal estimates; V: var-cov matrix of causal estimates
```
The `gwp` object contains three columns with names `SpleioP`, `JointExposuresP`, and `OutcomeP` which respectively represent P-values for horizontal pleiotropy testing with $S_\text{pleio}$ , a two degree of freedom joint chi-square test for the exposures, and a one degree of freedom chi-square test for the outcome. Where $\beta^Z_k$ generically denotes an effect size of association between SNP *k* and the trait *Z*, the P-values in `JointExposuresP` are from testing the null hypothesis $H_{0k}: \beta_k^\text{hbA1c}=\beta_k^\text{uric acid}=0$ and P-values in `OutcomeP` are from testing the null hypothesis $H_{0k}: \beta_k^\text{CAD}=0$. Where $(\theta^\text{HbA1c},\theta^\text{uric acid})$ denotes the causal effects of HbA1c and uric acid on CAD, P-values in `SpleioP` test the null hypothesis $H_{0k}: \beta_k^\text{CAD}-\beta_k^\text{HbA1c}\theta^\text{HbA1c}-\beta_k^\text{uric acid}\theta^\text{uric acid}=0$. See Lorincz-Comi & Yang et al. 2023 for more details.

