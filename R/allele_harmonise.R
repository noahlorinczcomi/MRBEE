allele_harmonise <- function(ref_panel, gwas_data) {

  # Make sure the reference panel has columns SNP, A1 and A2
  stopifnot(all(c("SNP", "A1", "A2") %in% colnames(ref_panel)))

  # Make sure the GWAS data has columns SNP, A1, A2, and Tstat
  stopifnot(all(c("SNP", "A1", "A2", "Zscore") %in% colnames(gwas_data)))

  # Merge the reference panel with the GWAS data
  merged_data <- merge(gwas_data, ref_panel, by = "SNP")

  # Remove any SNPs where A1 and A2 do not match between the reference panel and the GWAS data
  ind=which((toupper(merged_data$A1.x) == toupper(merged_data$A1.y) & toupper(merged_data$A2.x) == toupper(merged_data$A2.y))
            |(toupper(merged_data$A1.x) == toupper(merged_data$A2.y) & toupper(merged_data$A2.x) == toupper(merged_data$A1.y)))

  merged_data=merged_data[ind,]

  # Multiply Tstat by -1 if A1 and A2 are opposed between the reference panel and the GWAS data
  merged_data$Zscore <- ifelse(toupper(merged_data$A1.x) == toupper(merged_data$A2.y) & toupper(merged_data$A2.x) == toupper(merged_data$A1.y), -1 * merged_data$Zscore, merged_data$Zscore)

  colnames(merged_data)[which(names(merged_data) == "A1.y")] <- "A1"
  colnames(merged_data)[which(names(merged_data) == "A2.y")] <- "A2"
  colnames(merged_data)[which(names(merged_data) == "A1.x")] <- "A1.orginal"
  colnames(merged_data)[which(names(merged_data) == "A2.x")] <- "A2.orginal"

  return(merged_data)
}
