# MRBEE
This software accompanies the following papers:

- Lorincz-Comi, N., Yang, Y., Li, G., & Zhu, X. (2022). MRBEE: A novel bias-corrected multivariable Mendelian Randomization method. *bioRxiv*.
- Yang, Y., Lorincz-Comi, N., Zhu, X. (2022). Bias-corrected estimating equation of causal effect in multivariable Mendelian Randomization. *bioRxiv*.

Two pieces of software are provided in this repository:
- **MRBEE R package**
  - Install using `devtools::install_github("noahlorinczcomi/MRBEE")` or `remotes::install_github("noahlorinczcomi/MRBEE")` in R.
- **corrMatrix.py**
  - MRBEE subtracts from IVW two terms, both which are calculated from a correlation matrix **R**. If you have $p$ exposures and $q$ outcomes in MR, **R** will be of dimension $(p+q)\times (p+q)$ and the **MRBEE** software needs it.
  - **corrMatrix.py** is a command line tool to calculate **R** in a simple way (see below for example).
  - can be downloaded directly from repository.

# Examples
## Calculating **R** with **corrMatrix.py**
On your machine, download the **corrMatrix.py** file in a new directory named `/newdir`. In `/newdir`, move files containing GWAS summary statistics for all exposures and outcomes with which you intend to perform MR. 

## Performing MR with MRBEE
Here is an example of how to use MRBEE software for Mendelian Randomization with $m$ instrumental variables, $p$ exposures, and $q$ outcomes with fake generated data.

```
################################################### generating fake data
m=64;p=4;q=4; # number IVs, number exposures, number outcomes
n0=10000 # GWAS sample sizes of outcomes
n1=10000 # GWAS sample sizes of exposures
by=matrix(rnorm(m*q),nrow=m) # outcome GWAS estimates
bx=matrix(rnorm(m*p),nrow=m) # exposure GWAS estimates
sy=matrix(rchisq(m*q,n0-1)/n0,nrow=m) # outcome GWAS standard errors
sx=matrix(rchisq(m*p,n1-1)/n1,nrow=m) # exposure GWAS standard errors
mn=min(c(n0,n1)) 
R=rWishart(1,mn,diag(p+q))/mn # correlation matrix due to sample overlap and phenotypic correlation (calculated with corrMatrix.py or manually)
dim(R)=c(p+q,p+q)
cn1=paste0("x",1:p);cn2=paste0("y",1:q)
colnames(bx)=colnames(sx)=cn1 # assigning identical column names to `bx` and `sx` (required)
colnames(by)=colnames(sy)=cn2 # assigning identical column names to `by` and `sy` (required)
rownames(R)=colnames(R)=c(cn1,cn2) # assinging row, column names to R that are from column names of `bx` and `by` 
# (R can be arranged in any way, so long as its row and column names are correctly specified

################################################### performing analysis
pd=prepData(bx,by,sx,sy,R) # prepare data for MRBEE
fit=MRBEE(pd) # estimate causal effects using MRBEE 
Sp=Spleio(pd,fit$CausalEstimates,fit$VCovCausalEstimates) # Spleio statistics and P-values for horizontal pleiotropy for each IV 
```
