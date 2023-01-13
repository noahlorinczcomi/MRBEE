# MRBEE
This software accompanies the following papers:

- Lorincz-Comi, N., Yang, Y., Li, G., & Zhu, X. (2022). MRBEE: A novel bias-corrected multivariable Mendelian Randomization method. *bioRxiv*, doi.org/10.1101/2023.01.10.523480
- Yang, Y., Lorincz-Comi, N., Zhu, X. (2022). Bias-corrected estimating equation of causal effect in multivariable Mendelian Randomization. *arXiv*, https://doi.org/10.48550/arXiv.2301.05130

Please feel free to email Noah Lorincz-Comi (noahlorinczcomi@gmail.com, njl96@case.edu) with any questions.

Two pieces of software are provided in this repository:
- **MRBEE R package**
  - Install using `devtools::install_github("noahlorinczcomi/MRBEE")` or `remotes::install_github("noahlorinczcomi/MRBEE")` in R.
- **corrMatrix.py**
  - MRBEE subtracts from IVW two terms, both which are calculated from a correlation matrix **R**. If you have $p$ exposures and $q$ outcomes in MR, **R** will be of dimension $(p+q)\times (p+q)$ and the **MRBEE** software needs it.
  - **corrMatrix.py** is a command line tool to calculate **R** in a simple way (see below for example).
  - can be downloaded directly from repository.
    - ***If you do not want to use **corrMatrix.py**, you can do the following***:
      - Load all full exposure and outcome GWAS summary statistics into R
      - Standardize all association estimates (e.g., using the methods of Qi & Chatterjee, 2019 *Nature communications*)
      - Remove all SNPs with $P<\tau$ for any exposure or outcome (we recommend $\tau=0.05$)
      - Use the `cor()` function to calculate the correlations between all standardized exposure and outcome association estimates
      - The output of this code is the **R** that you need

# Examples
## Calculating **R** with **corrMatrix.py**
The **corrMatrix.py** program takes the following arguments:
- `-data`: (Required) comma-separated list of filepaths to all GWAS sumstats
- `-snp`: (Required) comma-separated list of SNP-identifying (i.e., rsID) column names for all GWAS sumstats (order must correspond to the order of the list passed to `-data`)
- `-beta`: (Required) comma-separated list of names of BETA columns in GWAS sumstats (order matters, always)
- `-se`: (Required) comma-separated list of names of SE columns in GWAS sumstats
- `-pt`: (Optional) P-value threshold. Only SNPs with all GWAS estimates with P>this threshold will be used to calculate **R**.
  - We recommend using `-pt 0.05`, which is the default value
- `-names`: (Required) Comma-separated list of row, column names to assign to the correlation matrix **R**
  - order corresponds to the order of data passed to `-data`
- `-out`: (Required) Filepath (extension optional, will be space-separated) of location to write out correlation matrix **R**
  - The correlation matrix **R** will also be printed on the screen

On your machine, download the **corrMatrix.py** file in a new directory named `/newdir`. Move all files containing GWAS summary statistics for all exposures and outcomes with which you intend to perform MR to `/newdir`. For the purpose of example, I have created 4 sets of simulated GWAS summary statistics, corresponding to two exposures and two outcomes, and put them in the `/newdir` directory. To demonstrate the flexibility of **corrMatrix.py**, the four GWAS summary statistic data sets have file extensions of ".txt", ".csv", ".txt.gz", and ".csv.gz".

Here is how you use **corrMatrix.py** with these simulated data:
```
cd /newdir
pwd
[1] gwasOutcomeset1.csv.gz  gwasOutcomeset2.txt  gwasExposureset1.txt.gz  gwasExposureset2.csv  corrMatrix.py
python corrMatrix.py \\
 -data gwasOutcomeset1.csv.gz,gwasOutcomeset2.txt,gwasExposureset1.txt.gz,gwasExposureset2.csv \\ 
 -snp rsID,rsID,rsID,rsID \\ 
 -beta betaOutcome1,betaOutcome2,betaExposure1,betaExposure2 \\ 
 -se seOutcome1,seOutcome2,seExposure1,seExposure2 
 -pt 0.05 \\
 -names y1,y2,x1,x2 \\
 -out R
```
which produces the following output
```
NOTE: -snp, -beta, -se, and -names flag declarations must be in the corresponding order of -data declarations
the program is running
          y1        y2        x1        x2
y1  1.000000  0.001672 -0.003739  0.000473
y2  0.001672  1.000000 -0.002589 -0.000470
x1 -0.003739 -0.002589  1.000000  0.005381
x2  0.000473 -0.000470  0.005381  1.000000
81038 SNPs used in correlation matrix estimation
```
and took approximately 10 seconds to run for 1,000,000 SNPs. Now, the correlation matrix **R** is saved in a space-delimited file named `R` in the `/newdir` directory.

## Performing MR with MRBEE
Here is an example of how to use MRBEE software for Mendelian Randomization with $m=100$ instrumental variables, $p=2$ exposures, and $q=2$ outcomes with fake generated data. We will use the matrix **R** calculated above using **corrMatrix.py** here.

```R
################# generating fake data #################
setwd("/newdir")
R=read.table("R", sep=" ",header=TRUE,row.names=1); R=as.matrix(R)
m=100;p=2;q=2; # number IVs, number exposures, number outcomes
n0=10000 # GWAS sample sizes of outcomes
n1=10000 # GWAS sample sizes of exposures
by=matrix(rnorm(m*q),nrow=m) # outcome GWAS estimates
bx=matrix(rnorm(m*p),nrow=m) # exposure GWAS estimates
sy=matrix(rchisq(m*q,n0-1)/n0,nrow=m) # outcome GWAS standard errors
sx=matrix(rchisq(m*p,n1-1)/n1,nrow=m) # exposure GWAS standard errors
cn1=paste0("x",1:p);cn2=paste0("y",1:q)
colnames(bx)=colnames(sx)=cn1 # assigning identical column names to `bx` and `sx` (required)
colnames(by)=colnames(sy)=cn2 # assigning identical column names to `by` and `sy` (required)
# NOTE: row and column names of R must contain all column names of `bx` and `by` 
# (R can be arranged in any way, so long as its row and column names are correctly specified)

################# performing analysis #################
library(MRBEE)
pd=prepData(bx,by,sx,sy,R) # prepare data for MRBEE
fit=MRBEE(pd) # estimate causal effects using MRBEE 
Sp=Spleio(pd,fit$CausalEstimates,fit$VCovCausalEstimates) # Spleio statistics and P-values for horizontal pleiotropy for each IV 
```
