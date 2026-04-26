# MRBEE

The MRBEE package is designed for conducting multivariable Mendelian Randomization (MVMR) analyses.

The full set of causal effect estimators in `MRBEE` is:
- `MRBEE.IMRP()`: MRBEE for MVMR
- `MRBEE.IMRP.Egger()`: MRBEE for MVMR with an intercept term
- `MRBEE.IMRP.UV()`: MRBEE for univariable (single-exposure)  MR

> [!NOTE]
> ```MRBEE.IMRP()``` removes weak instrument bias from GWAS estimation error from MR by using an unbiased estimating equation. MRBEE removes horizontal pleiotropy using a pleiotropy test applied iteratively during causal estimation.

## Installation

```R
devtools::install_github("noahlorinczcomi/MRBEE")
```

or

```R
remotes::install_github("noahlorinczcomi/MRBEE")
```

## Running (**TL;DR**)

```R
mrbee_result = MRBEE::MRBEE.IMRP(
  bX = <EFFECT SIZE EXPOSURES>,         # m x p matrix (m IVs, p exposures)
  by = <EFFECT SIZE OUTCOME>,           # m x 1 vector (m IVs)
  bXse = <SE EXPOSURES>,                # m x p matrix (m IVs, p exposures)
  byse = <SE OUTCOME>,                  # m x 1 vector (m IVs)
  Rxy = <ESTIMATION ERROR VCOV MATRIX>  # (p + 1) x (p + 1) matrix (p exposures, 1 outcome)
)

"""
NOTE:

1. `<ESTIMATION ERROR VCOV MATRIX>` can be calculated from (equivalently):
  - `MRBEE::errorcov()`
  - `ldscR::ldscR()`

2. MVMR IVs can be inferred by:
  a. Applying `MRBEE::Joint.test()` to `(<EFFECT SIZE EXPOSURES>, <SE EXPOSURES>)`
  b. LD clumping + significance thresholding (e.g., using `plink`)

3. `EFFECT SIZE <X>` / `SE <X>` can either be GWAS beta/standard error pairs, or Z-score/1 pairs

4. The final row/column of `Rxy` corresponds to the MR outcome; the rest to the MR exposures
"""
```

### Output
The output of `MRBEE::MRBEE.IMRP()` is a list with these named elements (classes):

- `theta` (`numeric`): causal effect estimates
- `covtheta` (`matrix`, `array`): variance-covariance matrix of causal effect estimates
- `delta` (`numeric`): residual; an estimate of horizontal pleiotropy

## Running (**worked example**)
We have pre-saved the following data to use in the example ([`example.R`](example.R)) that will follow:
- `http://tinyurl.com/nhdfwd8v`: Exposure and outcome GWAS summary statistics
- `http://tinyurl.com/8yt7bhvp`: Pre-identified instrumental variables to use in MVMR


### Data Preparation and Harmonization
First, download and prepare your GWAS summary data. Then, harmonize effect alleles between all exposure and outcome GWAS and the LD reference panel using the ```MRBEE::filter_align()``` function:

```R
# download example data (list of GWAS summary statistics)
data_url = "http://tinyurl.com/nhdfwd8v"
temp_file = tempfile()
download.file(data_url, temp_file, mode="wb")
gwaslist = readRDS(temp_file)
unlink(temp_file)

# perform MRBEE
library(MRBEE)
data("hapmap3") # LD reference panel
gwaslist = filter_align(
  gwas_data_list=gwaslist,
  ref_panel=hapmap3[,c("SNP", "A1", "A2")]
)
```

### Estimation Error Covariance Matrix
You can calculate the estimation error variance-covariance matrix using insignificant SNPs ([Zhu et al., 2015](https://doi.org/10.1016/j.ajhg.2014.11.011)) with `MRBEE::errorCov` like this:

```R
# matrix Z-statistics from all exposures and the outcome
ZMatrix = cbind(
  # exposure GWAS
  gwaslist$driving$Zscore,
  gwaslist$computer$Zscore,
  gwaslist$TV$Zscore,
  gwaslist$schooling$Zscore,
  # outcome GWAS
  gwaslist$myopia$Zscore
)

# estimation error variance-covariance matrix
Rxy = errorCov(ZMatrix=ZMatrix)
```

> [!CAUTION]
> The estimation error variance-covariance matrix `Rxy` is assumed by `MRBEE` to have the final row/column correspond to the **outcome**. For example, in MVMR with `p` exposures, `Rxy` is `(p + 1) x (p + 1)`, and `Rxy[p + 1, p + 1]` is the **outcome** estimation error variance. This makes `Rxy[1:p, 1:p]` the **exposure** estimation error variance-covariance matrix.

To alternatively use LD score regression ([Bulik-Sullivan et al., 2015](https://doi.org/10.1038/ng.3211)) from the `ldscR` package, run:

```R
devtools::install_github("harryyiheyang/ldscR")
library(ldscR)
data("EURLDSC")
Rxy = ldscR(GWAS_List=gwaslist, LDSC=EURLDSC)$ECovEst
```

> [!NOTE]
> Both ```ldscR``` and the insignificant GWAS effect estimation method ([Zhu et al., 2015](https://doi.org/10.1016/j.ajhg.2014.11.011)) can be used for estimating the covariance matrix. Slight numerical differences will exist between the two approaches, but these differences will not significantly impact the final results.

### Select MRBEE instrumental variables

Run:

```R
jointtest = Joint.test(
  bZ=ZMatrix[, -5], # remove MVMR outcome by index
  RZ=Rxy[-5,-5]     # remove MVMR outcome by index
)
jointtest$SNP = gwaslist$driving$SNP
```

> [!TIP]
> Instrumental variables (IVs) in multivariable MR (MVMR) can be selected from a chi-square joint test of association SNP-wise in the exposure set. This test can be more powerful than multiple exposure-specific tests, and can be applied using the `MRBEE::Joint.test()` function. **These test results are in general still subject to linkage disequilibrium (LD)**.

Now, `jointtest` contains joint test results for *all* SNPs. In the set of joint test-significant SNPs, there will be LD. A popular method used to select only independent SNPs from this set is "clumping and thresholding" (C+T) (e.g., using [`plink`](https://www.cog-genomics.org/plink/)).

To apply the C+T method implemented by `plink`, run:

```R
# PARAMS
ld_ref_panel = "data/1000G/1kg_phase3_EUR_only"
joint_test_results_file = "myopia/plinkfile/joint.txt"
c_plus_t_result_file = "myopia/plinkfile/joint"

# write SNP-P-value pairs to be read by PLINK cli
write.table(
  jointtest,
  joint_test_results_file,
  row.names=FALSE, quote=FALSE, sep="\t"
)

# run plink cli to perform C+T
system(
  paste(
    "plink --bfile",
    ld_ref_panel,
    "--clump",
    joint_test_results_file,
    "--clump-field P --clump-kb 500 --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --out",
    c_plus_t_result_file
  )
)

# read in C+T results (independent significant SNPs)
plink_ivs = data.table::fread(
  paste0(c_plus_t_result_file, ".clumped")
)$SNP
```

The pre-computed set of independent and joint test-significant SNPs (MVMR IVs) in our example can be loaded like this:

```R
plink_ivs = data.table::fread("http://tinyurl.com/8yt7bhvp", header = F)$V1
```

### Perform MVMR using `MRBEE`

> [!TIP]
> We recommend using Z-scores for MR analyses. If effect sizes are estimated using the same sample sizes for all SNPs and GWAS, MVMR estimates using Z-scores and effect sizes are equivalent. However, for sample sizes that vary across SNPs, Z-scores naturally assign sample size-proportional weight when estimating causal effects using MR. This resembles a second-stage reweighting, leading to more stable causal estimates.

To run `MRBEE` using our continued example with Z-scores, run:

```R
ZMatrix1 = ZMatrix[which(jointtest$SNP %in% plink_ivs), ]
fit_zscores = MRBEE.IMRP(
  by=ZMatrix1[, 5],                   # exposure Z-scores
  bX=ZMatrix1[, -5],                  # outcome Z-scores
  bXse=matrix(1, nrow(ZMatrix1), 4),  # exposure SE matrix (1s b/c using Z-scores)
  byse=rep(1, nrow(ZMatrix1)),        # outcome SE vector (1s b/c using Z-scores)
  Rxy=Rxy,                            # estimation error variance-covariance matrix
  var.est="ordinal"                   # residual variance estimate using delta method
)

print(fit_zscores$theta / sqrt(diag(fit_zscores$covtheta)))
```

### Perform MVMR using `MRBEE` with raw effect sizes and SEs

If you do not want to use Z-scores in MVMR, you can use raw GWAS effect size estimates and their SEs in `MRBEE` like this in our example:

```R
# SNP effect size estimates
BETA = cbind(
    # exposures
    gwaslist$driving$BETA,
    gwaslist$computer$BETA,
    gwaslist$TV$BETA,
    gwaslist$schooling$BETA,
    # outcome
    gwaslist$myopia$BETA)

# SNP standard errors
SE = cbind(
  # exposures
  gwaslist$driving$SE,
  gwaslist$computer$SE,
  gwaslist$TV$SE,
  gwaslist$schooling$SE,
  # outcome
  gwaslist$myopia$SE
)

# harmonise `BETA` like `Zmatrix` was harmonised above
# (the filter_align function only adjust the signs of Zscore)
BETA = abs(BETA) * sign(ZMatrix)

# filter `BETA` and `SE` to SNPs in the pre-identified IV set
BETA1 = BETA[which(jointtest$SNP %in% plink_ivs), ]
SE1 = SE[which(jointtest$SNP %in% plink_ivs), ]

# perform MRBEE
fit_betas = MRBEE.IMRP(
  by=BETA1[, 5],     # exposure effect sizes
  bX=BETA1[, -5],    # outcome effect sizes
  bXse=SE1[, -5],    # exposure SEs
  byse=SE1[, 5],     # outcome SEs
  Rxy=Rxy,           # estimation error variance-covariance matrix
  var.est="ordinal"  # residual variance estimate using the delta method
)

print(fit_betas$theta / sqrt(diag(fit_betas$covtheta)))
```

## Metadata
MRBEE publication
- Lorincz-Comi, N., Yang, Y., Li, G., & Zhu, X. (2024). MRBEE: a bias-corrected multivariable Mendelian randomization method. Human Genetics and Genomics Advances, 5(3). [https://doi.org/10.1016/j.xhgg.2024.100290](https://doi.org/10.1016/j.xhgg.2024.100290)

Maintainers
- Yihe Yang (`yxy1234@case.edu`, ORCID: `0000-0001-6563-3579`)
- Noah Lorincz-Comi (`njl96@case.edu`, ORCID: `0000-0002-0517-2499`)
