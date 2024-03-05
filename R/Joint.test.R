#' Perform Joint Test on GWAS Summary Data
#'
#' This function conducts a joint test on GWAS summary data, calculating Chi-squared values
#' and corresponding P-values based on allele-aligned Z-scores and the correlation matrix
#' of estimation errors for these Z-scores.
#'
#' @param bZ A matrix of allele-aligned Z-scores from GWAS summary data, where each row
#'   represents a different genetic variant and each column represents a different trait
#'   or study. The matrix dimensions are n x p, with n being the number of variants and
#'   p being the number of traits or studies.
#' @param RZ The correlation matrix of estimation errors for the Z-scores in `bZ`. This
#'   matrix should be p x p, where p is the number of traits or studies.
#'
#' @return A data frame with two columns: `Chi2`, containing the Chi-squared values for
#'   each genetic variant, and `P`, containing the corresponding P-values.
#'
#' @examples
#' # Example usage:
#' # Assuming bZ is your matrix of allele-aligned Z-scores and RZ is the correlation matrix
#' result <- Jointtest(bZ, RZ)
#'
#' @export
Joint.test <- function(bZ, RZ) {
  n <- nrow(bZ)
  p <- ncol(bZ)
  z <- pv <- bZ[1,]
  ThetaZ <- solve(RZ)
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in 1:n) {
    a <- c(bZ[i,])
    b <- c(ThetaZ %*% a)
    chi <- sum(a * b)
    z[i] <- chi
    setTxtProgressBar(pb, i)
    pv[i] <- stats::pchisq(chi, p, lower.tail = F)
  }
  close(pb)
  return(data.frame(Chi2 = z, P = pv))
}
