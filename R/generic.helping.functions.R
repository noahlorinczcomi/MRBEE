#' helping functions
 
array.to.list=function(.array) lapply(seq(dim(.array)[3]), function(x) .array[ , , x])

ssum=function(.list) {
  if(is.array(.list)) .list=array.to.list(.list)
  return(Reduce("+", .list))
}

subset.all=function(A, B, SigmaUU, SigmaUV, SigmaVV, Outliers) {
	if(is.null(Outliers) | length(Outliers) == 0) {
		A=as.matrix(A); B=as.matrix(B)
		return(list(A=A, B=B, SigmaUU=SigmaUU, SigmaUV=SigmaUV, SigmaVV=SigmaVV))
	} else {
		A=as.matrix(A); A=A[-Outliers,]; B=as.matrix(B); B=B[-Outliers,]
		SigmaUU=SigmaUU[,,-Outliers, drop=FALSE]
		SigmaVV=SigmaVV[,,-Outliers,drop=FALSE]
		SigmaUV=SigmaUV[,,-Outliers,drop=FALSE]
		return(list(A=A, B=B, SigmaUU=SigmaUU, SigmaUV=SigmaUV, SigmaVV=SigmaVV))
	}
}