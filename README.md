# MRBEE
This software accompanies the following paper:

- Lorincz-Comi, N., Yang, Y., Li, G., & Zhu, X. (2023). MRBEE: A novel bias-corrected multivariable Mendelian Randomization method. *bioRxiv*, https://doi.org/10.1101/2023.01.10.523480

Please feel free to email Noah Lorincz-Comi (noahlorinczcomi@gmail.com, njl96@case.edu) with any questions.

## TL;DR
```R
# bx: mxp matrix of standardized IV associations with exposures
# bxse: mxp matrix of standardized IV SEs for each exposures
# by: mx1 vector of standardized IV associations with outcome
# byse: mx1 vector of standardized IV SEs for the outcome
# R: (p+1)x(p+1) matrix of correlations between measurement errors for the outcome (first/top left position) and each exposure
# Ncor: number of nonsignificant SNPs used to calculate `R` (if known to be large, eg >10,000, you can set this to be any number >10,000)
bT=list(R=R,Ncor=Ncor,EstHarm=cbind(by,bx),SEHarm=cbind(byse,bxse))
pD=prepData(bT)
fit=MRBEE.IMRP(pD) # stores causal estimates and some model characteristics
res=data.frame(Est=fit$CausalEstimates,SE=sqrt(diag(fit$VCovCausalEstimates))); res$P=1-pchisq((res$Est/res$SE)^2,1)
res # simplified results
Sp=Spleio(pd,fit$CausalEstimates,fit$VCovCausalEstimates) # Spleio statistics and P-values for horizontal pleiotropy for each IV 
```

<!---
Two pieces of software are provided in this repository:
- **MRBEE R package**
  - Install using `devtools::install_github("noahlorinczcomi/MRBEE")` or `remotes::install_github("noahlorinczcomi/MRBEE")` in R.
- `corrMatrix.py`
  - MRBEE subtracts from IVW two terms, both of which are calculated from a correlation matrix **R**. If you have $p$ exposures and $q$ outcomes in MR, **R** will be of dimension $(p+q)\times (p+q)$ and the **MRBEE** software needs it.
  - `corrMatrix.py` is a command line tool to calculate **R** in a simple way (see below for example).
  - ***If you do not want to use*** `corrMatrix.py`, ***you can do the following***:
    - Load all full exposure and outcome GWAS summary statistics into R
    - Standardize all association estimates (e.g., using the methods of Qi & Chatterjee, 2019 *Nature communications*)
    - Remove all SNPs with $P<\tau$ for any exposure or outcome (we recommend $\tau=0.05$)
    - Use the `cor()` function to calculate the correlations between all standardized exposure and outcome association estimates
    - The output of this code is the **R** that you need

# Example
## Calculating **R** with `corrMatrix.py`
The `corrMatrix.py` program takes the following arguments:
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

On your machine, download the `corrMatrix.py` file in a new directory named `/newdir`. Move all files containing GWAS summary statistics for all exposures and outcomes with which you intend to perform MR to `/newdir`. For the purpose of example, I have created 4 sets of simulated GWAS summary statistics, corresponding to two exposures and two outcomes, and put them in the `/newdir` directory. To demonstrate the flexibility of `corrMatrix.py`, the four GWAS summary statistic data sets have file extensions of ".txt", ".csv", ".txt.gz", and ".csv.gz".

Here is how you use `corrMatrix.py` with these simulated data:
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
and took approximately 10 seconds to run for 1,000,000 SNPs. Now, the correlation matrix $\mathbf{R}$ is saved in a space-delimited file named `R` in the `/newdir` directory.
--->

## Performing MR with MRBEE
MRBEE requires the matrix $\mathbf{R}$ of correlations between measurement errors for the outcome phenotype and all exposure phenotypes. In the code below, this is the matrix `R`. You can calculate this matrix using insignificant GWAS summary statistics or LD score regression (see the [MRBEE paper](https://doi.org/10.1101/2023.01.10.523480 ), [MVMR-cML paper](https://doi.org/10.1016/j.ajhg.2023.02.014), or [CPASSOC paper](https://doi.org/10.1016/j.ajhg.2014.11.011)). 

In the tutorial below, we will show you how to calculate $\mathbf{R}$ using the R language. You can also use our `corrMatrix.py` tool to calculate $\mathbf{R}$. See its tutorial [here](https://github.com/noahlorinczcomi/MRBEE/blob/main/corrMatrix_Tutorial.ipynb).

Below is an example of how to use MRBEE software for multivariable Mendelian Randomization. In this example, we use publicly available GWAS data for CAD and 9 cardiometabolic exposures.

### Downloading data
If your data is not already downloaded and ready to be analyzed, you can use data that we provide to complete this tutorial. First, download these data using the following commands in a Linux terminal:
```unix
wget http://hal.case.edu/~njl96/cardioData.txt.gz
wget http://hal.case.edu/~njl96/ivrsids.txt
```
Note that even if you are unable to download these data, you may still be able to read the tutorial to see how to use **MRBEE** with your own data.
### Preparing data
We can now load the data in R and look at its column names:
```R
IVrsIDs=data.table::fread("ivrsids.txt")
cardioData=data.table::fread("cardioData.txt.gz")
[1] "rsID" "CHR_BP" "BETA_CAD" "BETA_HDL" "BETA_LDL" "BETA_TG" "BETA_SBP"
[2] "BETA_UA" "BETA_HBA1C" "BETA_BMI" "BETA_HEIGHT" "BETA_HG" "SE_CAD"
[3] "SE_HDL" "SE_LDL" "SE_TG" "SE_SBP" "SE_UA" "SE_HBA1C" "SE_BMI"
[4] "SE_HEIGHT" "SE_HG" "EFFECT_ALLELE_CAD" "EFFECT_ALLELE_HDL"
[5] "EFFECT_ALLELE_LDL" "EFFECT_ALLELE_TG" "EFFECT_ALLELE_SBP"
[6] "EFFECT_ALLELE_UA" "EFFECT_ALLELE_HBA1C" "EFFECT_ALLELE_BMI"
[7] "EFFECT_ALLELE_HEIGHT" "EFFECT_ALLELE_HG" "NONEFFECT_ALLELE_CAD"
[8] "NONEFFECT_ALLELE_HDL" "NONEFFECT_ALLELE_LDL" "NONEFFECT_ALLELE_TG"
[9] "NONEFFECT_ALLELE_SBP" "NONEFFECT_ALLELE_UA" "NONEFFECT_ALLELE_HBA1C"
[10] "NONEFFECT_ALLELE_BMI" "NONEFFECT_ALLELE_HEIGHT" "NONEFFECT_ALLELE_HG"
[11] "MAF_CAD"
```
For simplicity, we create a vector of exposure names using:
```R
exposures=c("HDL","LDL","TG","SBP","UA","HBA1C","BMI","HEIGHT","HG")
```
Next, we need to specify the column names of beta, SE, and effect alleles for the outcome and each exposure like this:
```R
ests=c("BETA_CAD",paste0("BETA_",exposures"))
ses=c("SE_CAD",paste0("SE_",exposures))
alleles=c("EFFECT_ALLELE_CAD",paste0("EFFECT_ALLELE_", exposures))
```
It is required that the $k$th index position of `ests` corresponds to the same phenotype as the $k$th index positions of `ses` and `alleles`. A column naming convention such as the one above of the format `BETA_<x>`, `SE_<x>`, and `EFFECT_ALLELE_<x>` ensures that this requirement is met without much work, although the user is free to choose any naming convention so long as `ests`, `ses`, and `alleles` are specified with the correct orderings of their elements. We consolidate this information by putting each into a single list:
```R
dNames=list(est=ests,se=ses,allele=alleles)
```
### Calculating bias-correction terms
We calculate the bias-correction terms that MRBEE uses to correct for bias in multivariable IVW:
```R
bT=biasTerms(cardioData,dNames)
names(bT)
[1] "R" "Ncor" "EstHarm" "SEHarm"
```
These elements contain the following information:
1. `R`: $(p+1)\times(p_1)$ correlation matrix between measurement errors for the outcome and all exposures
2. `Ncor`: number of nonsignificant SNPs that were used to calculate `R`
3. `EstHarm`: Allele-harmonised association estimates for all SNPs, which the user should have already put in standardized scale
4. `SEHarm`: Standard errors corresponding to values in `EstHarm`, which the user should have already put in standardized scale

In the example so far, we have put column name values corresponding to the outcome at the front of each vector in which we are required to specify them (e.g., `ests`, `se`, `alleles`. This is a good practice and MRBEE assumes this by default. If the index position of values corresponding to the outcome are in a different position, set `oi=x` where `x` is the correct index position.

```R
oi=1 # this column of `EstHarm` and `SEHarm` correspond to the outcome if `ests`[1] corresponded to the outcome
IVInds=which(cardioData$rsID %in% IVrsIDs$SNP) # IVIinds now contains index positions of IVs to use
pD=prepData(bT,IVInds,outcome_index=oi)
names(pD)
[1] "betaX" "betaY" "UU" "UV" "VV"
```
Now the values in `pD` are subsetted to only those corresponding to the IVs. Now we can perform MVMR with MRBEE:
```R
fit=MRBEE.IMRP(pD)
names(fit)
[1] "CausalEstimates" "VCovCausalEstimates" "CausalEstimatePS"
[4] "InitialCausalEstimates" "VcovInitialCausalEstimates" "PleiotropyIndices"
[7] "nIterations" "Qstart" "Qend" "ThetaDifferences"
```
