#' Estimate Error Covariance Matrix Using GWAS Insignificant Effects
#'
#' This function estimates the error covariance matrix by subsampling a proportion of insignificant GWAS effects and calculating their correlation coefficients.
#'
#' @param ZMatrix A matrix of Z-scores for exposure and outcome, with the outcome GWAS in the last column.
#' @param Zscore.cutoff The cutoff for significance. Defaults to 2.
#' @param subsampling.ratio The proportion of effects to subsample for each iteration. Defaults to 0.1.
#' @param subsampling.time The number of subsampling iterations. Defaults to 1000.
#' @return A matrix representing the estimated error covariance.
#' @export
#'
errorCov=function(ZMatrix,Zscore.cutoff=2,subsampling.ratio=0.1,subsampling.time=1000){
z=apply(abs(ZMatrix),1,max)
insignind=which(z<=Zscore.cutoff)
p=ncol(ZMatrix)
n=length(insignind)
RList=array(0,c(subsampling.time,p,p))
pb <- utils::txtProgressBar(min = 0, max = subsampling.time, style = 3)
for(i in 1:subsampling.time){
ind=sample(insignind,subsampling.ratio*n)
Z=ZMatrix[ind,]
RList[i,,]=stats::cor(Z)
utils::setTxtProgressBar(pb, i)
}
close(pb)
R=RList[1,,]
for(i in 1:p){
for(j in i:p){
R[i,j]=R[j,i]=stats::median(RList[,i,j])
}
}
return(R)
}
