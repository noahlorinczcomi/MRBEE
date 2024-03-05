IVweight=function(byse,bXse,Rxy){
  bZse=cbind(bXse,byse)
  p=dim(bZse)[2]
  n=dim(bZse)[1]
  RxyList=array(0,c(n,p,p))
  for(i in 1:n){
    s=bZse[i,]
    RxyList[i,,]=t(t(Rxy)*s)*s
  }
  return(RxyList)
}

imrpdetect=function(x,theta,RxyList,indvalid,var.est="robust",FDR=T,adjust.method="Sidak"){
  p=length(theta)
  if(var.est=="robust"){
    varx=stats::mad(x[indvalid])^2
  }
  if(var.est=="variance"){varx=stats::var(x[indvalid])}
  if(var.est=="ordinal"){
    varx=x*0
    for(i in 1:length(x)){
      varx[i]=c(RxyList[i,p+1,p+1]+t(theta)%*%RxyList[i,1:p,1:p]%*%theta-2*sum(theta*RxyList[i,p+1,1:p]))
    }
  }
  pv=stats::pchisq(x^2/varx,1,lower.tail=F)
  if(FDR==T){
    pv=p.fdr(pvalues=pv,adjust.method=adjust.method)$fdrs
  }
  return(as.vector(pv))
}

validadj <- function(vector1, vector2, tau) {
  diff <- length(vector2) / length(vector1)
  if (diff < tau) {
    missing_indices <- setdiff(1:length(vector1), vector2)
    sorted_missing_indices <- missing_indices[order(vector1[missing_indices])]
    num_to_add <- ceiling(tau * length(vector1)) - length(vector2)
    vector2 <- c(vector2, sorted_missing_indices[1:num_to_add])
  }
  return(vector2)
}

biasterm=function(RxyList,indvalid){
  X=RxyList[1,,]*0
  for(i in indvalid){
    X=X+RxyList[i,,]
  }
  return(X)
}
