test_that("allele_harmonise keeps Z-score when alleles match reference", {
  ref <- data.frame(SNP = "rs1", A1 = "A", A2 = "G", stringsAsFactors = FALSE)
  gwas <- data.frame(SNP = "rs1", A1 = "A", A2 = "G", Zscore = 1.5, stringsAsFactors = FALSE)
  res <- MRBEE:::allele_harmonise(ref, gwas)
  expect_equal(res$Zscore, 1.5)
})

test_that("allele_harmonise flips Z-score when alleles are reversed", {
  ref <- data.frame(SNP = "rs1", A1 = "A", A2 = "G", stringsAsFactors = FALSE)
  gwas <- data.frame(SNP = "rs1", A1 = "G", A2 = "A", Zscore = 1.5, stringsAsFactors = FALSE)
  res <- MRBEE:::allele_harmonise(ref, gwas)
  expect_equal(res$Zscore, -1.5)
})

test_that("allele_harmonise removes SNPs with incompatible alleles", {
  ref <- data.frame(SNP = "rs1", A1 = "A", A2 = "G", stringsAsFactors = FALSE)
  gwas <- data.frame(SNP = "rs1", A1 = "C", A2 = "T", Zscore = 1.5, stringsAsFactors = FALSE)
  res <- MRBEE:::allele_harmonise(ref, gwas)
  expect_equal(nrow(res), 0L)
})

test_that("allele_harmonise is case-insensitive", {
  ref <- data.frame(SNP = "rs1", A1 = "a", A2 = "g", stringsAsFactors = FALSE)
  gwas <- data.frame(SNP = "rs1", A1 = "G", A2 = "A", Zscore = 2.0, stringsAsFactors = FALSE)
  res <- MRBEE:::allele_harmonise(ref, gwas)
  expect_equal(res$Zscore, -2.0)
})

test_that("allele_harmonise errors when reference panel is missing required columns", {
  ref <- data.frame(SNP = "rs1", A1 = "A", stringsAsFactors = FALSE)
  gwas <- data.frame(SNP = "rs1", A1 = "A", A2 = "G", Zscore = 1.0, stringsAsFactors = FALSE)
  expect_error(MRBEE:::allele_harmonise(ref, gwas))
})

test_that("allele_harmonise errors when GWAS data is missing required columns", {
  ref <- data.frame(SNP = "rs1", A1 = "A", A2 = "G", stringsAsFactors = FALSE)
  gwas <- data.frame(SNP = "rs1", A1 = "A", A2 = "G", stringsAsFactors = FALSE)
  expect_error(MRBEE:::allele_harmonise(ref, gwas))
})
