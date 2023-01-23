#' A helping function to perform a global test of unbalanced horizontal pleiotropy
#'
#' @param C a contrast matrix created by the iter.estimator() function (called from within the MRBEE() function) to extract intercept-relevant terms from a matrix
#' @param Est MRBEE causal estimates with first row corresponding to intercept terms
#' @param Var Variance-covariance matrix of MRBEE causal estimates
#' @keywords contrast testing
#' @export 
#' @examples
#' fQ()

fQ=function(C, Est, Var) t(C %*% c(Est)) %*% solve(C %*% Var %*% t(C)) %*% (C %*% c(Est))