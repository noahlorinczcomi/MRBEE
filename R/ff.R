#' A helping function
#'
#' This function finds the index locations of phenotypes in a correlation matrix.
#' @param R any square matrix where row names match column names.
#' @param x matrix with column names to locate in R.
#' @keywords helping
#' @export
#' @examples ind=sapply(colnames(X), function(h) ff(R,h))
#' ff()

ff=function(R,x) sapply(x,function(h) which(x[h]==colnames(R)))
