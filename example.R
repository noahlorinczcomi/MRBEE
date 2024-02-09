data_url <- "http://tinyurl.com/nhdfwd8v"
temp_file <- tempfile()
download.file(data_url, temp_file, mode="wb")
gwaslist=readRDS(temp_file)
unlink(temp_file)
library(ldscR)
library(devtools)
install()
library(MRBEE.IMRP)
data("hapmap3")
data("EURLDSC")

############################ Merge the data and estimate the correlation matrix of estimation error using LDSC

gwaslist=merge_intersect(gwas_data_list=gwaslist,ref_panel=hapmap3[,c("SNP","A1","A2")])
fitldsc=ldscR(GWAS_List=gwaslist,LDSC=EURLDSC)
ZMatrix=cbind(gwaslist$driving$Zscore,gwaslist$computer$Zscore,gwaslist$TV$Zscore,gwaslist$schooling$Zscore,gwaslist$myopia$Zscore)
jointtest=Joint.test(bZ=ZMatrix[,-5],RZ=fitldsc$ECovEst[-5,-5])
jointtest$SNP=gwaslist$driving$SNP

########################### Perform C+T using PLINK ########################################

#write.table(jointtest,"myopia/plinkfile/joint.txt",row.names=F,quote=F,sep="\t")
#setwd("~/Plink")
#system("./plink --bfile data/1000G/1kg_phase3_EUR_only --clump myopia/plinkfile/joint.txt --clump-field P  --clump-kb 500 --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --out myopia/plinkfile/joint")
#plink=fread("~/myopia/plinkfile/joint.clumped")
plink=data.table::fread("http://tinyurl.com/4kbrzk32", header = TRUE)

########################### Using Zscore (We suggest) ######################################

ZMatrix1=ZMatrix[which(jointtest$SNP%in%plink$SNP),]
fit=MRBEE.IMRP(by=ZMatrix1[,5],bX=ZMatrix1[,-5],byse=rep(1,771),bXse=matrix(1,nrow(ZMatrix1),4),Rxy=fitldsc$ECovEst,var.est="ordinal")
print(fit$theta/sqrt(diag(fit$covtheta)))

############################## Using Effect Size and Standard Error ##################################

BETA=cbind(gwaslist$driving$BETA,gwaslist$computer$BETA,gwaslist$TV$BETA,gwaslist$schooling$BETA,gwaslist$myopia$BETA)
BETA=abs(BETA)*sign(ZMatrix) # The merge_intersect function only adjust the signs of Zscore
SE=cbind(gwaslist$driving$SE,gwaslist$computer$SE,gwaslist$TV$SE,gwaslist$schooling$SE,gwaslist$myopia$SE)
BETA1=BETA[which(jointtest$SNP%in%plink$SNP),]
SE1=SE[which(jointtest$SNP%in%plink$SNP),]
fit1=MRBEE.IMRP(by=BETA1[,5],bX=BETA1[,-5],byse=SE1[,5],bXse=SE1[,-5],Rxy=fitldsc$ECovEst,var.est="ordinal")
print(fit1$theta/sqrt(diag(fit1$covtheta)))
