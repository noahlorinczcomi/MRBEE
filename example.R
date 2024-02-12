data_url <- "http://tinyurl.com/nhdfwd8v"
temp_file <- tempfile()
download.file(data_url, temp_file, mode="wb")
gwaslist=readRDS(temp_file)
unlink(temp_file)
library(MRBEE)
data("hapmap3")

############################ Merge the data and estimate the correlation matrix of estimation error using LDSC
gwaslist=filter_align(gwas_data_list=gwaslist,ref_panel=hapmap3[,c("SNP","A1","A2")])
ZMatrix=cbind(gwaslist$driving$Zscore,gwaslist$computer$Zscore,gwaslist$TV$Zscore,gwaslist$schooling$Zscore,gwaslist$myopia$Zscore)
Rxy=errorCov(ZMatrix=ZMatrix)
# One can also use the LDSC implemented by our R package ldscR to estimate the
#library(ldscR)
#data("EURLDSC")
#fitldsc=ldscR(GWAS_List=gwaslist,LDSC=EURLDSC)
#Rxy=fitldsc$ECovEst
########################### Perform C+T using PLINK ########################################
jointtest=Joint.test(bZ=ZMatrix[,-5],RZ=Rxy[-5,-5])
jointtest$SNP=gwaslist$driving$SNP
#write.table(jointtest,"myopia/plinkfile/joint.txt",row.names=F,quote=F,sep="\t")
#setwd("~/Plink")
#system("./plink --bfile data/1000G/1kg_phase3_EUR_only --clump myopia/plinkfile/joint.txt --clump-field P  --clump-kb 500 --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --out myopia/plinkfile/joint")
#plink=fread("~/myopia/plinkfile/joint.clumped")[,1]
plink=data.table::fread("http://tinyurl.com/8yt7bhvp", header = F)

########################### Using Zscore (We suggest) ######################################
ZMatrix1=ZMatrix[which(jointtest$SNP%in%plink$V1),]
fit=MRBEE.IMRP(by=ZMatrix1[,5],bX=ZMatrix1[,-5],byse=rep(1,nrow(ZMatrix1)),bXse=matrix(1,nrow(ZMatrix1),4),Rxy=Rxy,var.est="ordinal")
print(fit$theta/sqrt(diag(fit$covtheta)))

############################## Using Effect Size and Standard Error ##################################

BETA=cbind(gwaslist$driving$BETA,gwaslist$computer$BETA,gwaslist$TV$BETA,gwaslist$schooling$BETA,gwaslist$myopia$BETA)
BETA=abs(BETA)*sign(ZMatrix) # The merge_intersect function only adjust the signs of Zscore
SE=cbind(gwaslist$driving$SE,gwaslist$computer$SE,gwaslist$TV$SE,gwaslist$schooling$SE,gwaslist$myopia$SE)
BETA1=BETA[which(jointtest$SNP%in%plink$V1),]
SE1=SE[which(jointtest$SNP%in%plink$V1),]
fit1=MRBEE.IMRP(by=BETA1[,5],bX=BETA1[,-5],byse=SE1[,5],bXse=SE1[,-5],Rxy=Rxy,var.est="ordinal")
print(fit1$theta/sqrt(diag(fit1$covtheta)))
