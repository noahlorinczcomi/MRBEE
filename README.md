# MRBEE.IMRP
The MRBEE.IMRP package is designed for conducting multivariable Mendelian Randomization (MVMR) analyses. MRBEE.IMRP removes weak instrument bias, which is caused by the estimation error of exposures and outcome GWAS, by using an unbiased estimatin function. On the other hand, MRBEE iteratively detected and remove the horiztonal pleiotropy using a pleiotropy test, making it robust to horizontal pleiotropy.

## Installation
You can install the MRBEE.IMRP package directly from GitHub using the following command:
```R
# Install the devtools package if you haven't already
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install MRBEE.IMRP from GitHub
devtools::install_github("harryyiheyang/MRBEE.IMRP")
```
Additionally, this package utilizes the ldscR package for allele harmonization and estimation of the covariance matrix of estimation errors. You can install ldscR from GitHub as well:
```R
# Install ldscR from GitHub
devtools::install_github("harryyiheyang/ldscR")
```

## Usage
Here's an example workflow using MRBEE.IMRP and ldscR:

### Data Preparation and Harmonization
First, download and prepare your GWAS summary data. Then, harmonize alleles using the merge_intersect() function from the ldscR package:
```R
library(ldscR)
library(MRBEE.IMRP)

# Example data download
data_url <- "http://tinyurl.com/nhdfwd8v"
temp_file <- tempfile()
download.file(data_url, temp_file, mode="wb")
gwaslist <- readRDS(temp_file)
unlink(temp_file)

# Harmonize data
data("hapmap3")
gwaslist <- merge_intersect(gwas_data_list=gwaslist, ref_panel=hapmap3[,c("SNP","A1","A2")])
```
### Estimation Error Covariance Matrix
```R
fitldsc <- ldscR(GWAS_List=gwaslist, LDSC=EURLDSC)
```
While both ldscR and the insignificant GWAS effect estimation method can be used for estimating the covariance matrix, slight differences exist between the two approaches. However, these differences do not significantly impact the final results.

### Perform Joint test to detect significant IVs
```R
ZMatrix=cbind(gwaslist$driving$Zscore,gwaslist$computer$Zscore,gwaslist$TV$Zscore,gwaslist$schooling$Zscore,gwaslist$myopia$Zscore)
jointtest=Joint.test(bZ=ZMatrix[,-5],RZ=fitldsc$ECovEst[-5,-5])
jointtest$SNP=gwaslist$driving$SNP

########################### Perform C+T using PLINK ########################################
#write.table(jointtest,"myopia/plinkfile/joint.txt",row.names=F,quote=F,sep="\t")
#setwd("~/Plink")
#system("./plink --bfile data/1000G/1kg_phase3_EUR_only --clump myopia/plinkfile/joint.txt --clump-field P  --clump-kb 500 --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --out myopia/plinkfile/joint")
#plink=fread("~/myopia/plinkfile/joint.clumped")
plink=data.table::fread("http://tinyurl.com/4kbrzk32", header = TRUE)
```
### Perform MVMR using Zscore
```R
ZMatrix1=ZMatrix[which(jointtest$SNP%in%plink$SNP),]
fit=MRBEE.IMRP(by=ZMatrix1[,5],bX=ZMatrix1[,-5],byse=rep(1,771),bXse=matrix(1,nrow(ZMatrix1),4),Rxy=fitldsc$ECovEst,var.est="ordinal")
print(fit$theta/sqrt(diag(fit$covtheta)))
```

### Perform MVMR using Effect Size and Standard Error
```R
BETA=cbind(gwaslist$driving$BETA,gwaslist$computer$BETA,gwaslist$TV$BETA,gwaslist$schooling$BETA,gwaslist$myopia$BETA)
BETA=abs(BETA)*sign(ZMatrix) # The merge_intersect function only adjust the signs of Zscore
SE=cbind(gwaslist$driving$SE,gwaslist$computer$SE,gwaslist$TV$SE,gwaslist$schooling$SE,gwaslist$myopia$SE)
BETA1=BETA[which(jointtest$SNP%in%plink$SNP),]
SE1=SE[which(jointtest$SNP%in%plink$SNP),]
fit1=MRBEE.IMRP(by=BETA1[,5],bX=BETA1[,-5],byse=SE1[,5],bXse=SE1[,-5],Rxy=fitldsc$ECovEst,var.est="ordinal")
print(fit1$theta/sqrt(diag(fit1$covtheta)))
```

## Maintainers
This package is maintained by:

Yihe Yang
Email: yxy1234@case.edu
ORCID: 0000-0001-6563-3579
