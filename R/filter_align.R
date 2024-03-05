#' Filter and Align GWAS Data to a Reference Panel
#'
#' The \code{filter_align} function processes a list of GWAS summary statistics data frames, harmonizes alleles according to a reference panel, removes duplicates, and aligns data to common SNPs. It's used to prepare data for further analysis such as LDSC.
#'
#' @param gwas_data_list A list of data.frames where each data.frame contains GWAS summary statistics for a trait. Each data.frame should include columns for SNP identifiers, Z-scores of effect size estimates, sample sizes (N), effect allele (A1), and reference allele (A2).
#' @param ref_panel A data.frame containing the reference panel data. It must include columns for SNP, A1, and A2.
#' @param allele_match Logical. Whether to match alleles. Default \code{TRUE}.
#'
#' @return A list of data.frames, each corresponding to an input GWAS summary statistics data frame, but filtered, harmonized, and aligned to the common SNPs found across all data frames.
#'
#' @examples
#' \dontrun{
#' # Assuming GWAS_List and ref_panel are already defined:
#' GWAS_List <- filter_align(GWAS_List, ref_panel)
#' }
#' @details The function performs several key steps: adjusting alleles according to a reference panel, removing duplicate SNPs, and aligning all GWAS data frames to a set of common SNPs. This is often a necessary preprocessing step before performing genetic correlation and heritability analyses.
#' @importFrom data.table setDT
#' @export

filter_align <- function(gwas_data_list, ref_panel, allele_match=TRUE) {

  cat("Adjusting effect allele according to reference panel...\n")
  p <- length(gwas_data_list)
  if(allele_match==T){
    for (i in 1:p) {
      A <- gwas_data_list[[i]]
      A <- setDT(A)
      A <- allele_harmonise(ref_panel = ref_panel, gwas_data = A)
      gwas_data_list[[i]] <- A
    }
  }
  gwas_data_list <- lapply(gwas_data_list, function(df) {
    df <- df[!duplicated(df$SNP), ]
    return(df)
  })

  cat("Finding common SNPs...\n")
  snp_sets <- lapply(gwas_data_list, function(df) {
    return(as.character(df$SNP))
  })

  common_snps <- Reduce(intersect, snp_sets)

  cat("Aligning data to common SNPs and ordering...\n")
  gwas_data_common_aligned <- lapply(gwas_data_list, function(df) {
    df_common <- df[df$SNP %in% common_snps, ]
    df_common <- df_common[order(df_common$SNP), ]
    rownames(df_common) <- 1:nrow(df_common)
    return(df_common)
  })

  cat("Filtering complete.\n")
  return(gwas_data_common_aligned)
}
