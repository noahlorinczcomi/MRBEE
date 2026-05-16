suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(data.table)
  library(MRBEE)
})

option_list <- list(
  make_option(c("-c", "--config"), type = "character", metavar = "FILE",
              help = "Path to analysis config YAML (see analyses/template.yaml)"),
  make_option(c("-o", "--out"), type = "character", default = "mrbee_results.csv",
              metavar = "FILE", help = "Output CSV path [default: %default]"),
  make_option("--var-est", type = "character", default = "variance",
              metavar = "METHOD",
              help = "Variance estimator: variance, robust, or ordinal [default: %default]"),
  make_option("--pv-thres", type = "double", default = 0.05,
              metavar = "NUM", help = "Pleiotropy p-value threshold [default: %default]"),
  make_option("--no-fdr", action = "store_true", default = FALSE,
              help = "Disable FDR correction for pleiotropy test")
)

parser <- OptionParser(
  usage = "%prog -c config.yaml [-o results.csv]",
  option_list = option_list
)
opt <- parse_args(parser)

if (is.null(opt$config)) {
  print_help(parser)
  stop("--config is required", call. = FALSE)
}

# ── 1. Config ────────────────────────────────────────────────────────────────
cat("[1/6] Reading config:", opt$config, "\n")
cfg <- yaml::read_yaml(opt$config)

if (is.null(cfg$ld) || !nzchar(cfg$ld)) {
  stop("'ld' must be set in the config (path to PLINK bfile prefix, without extension)",
       call. = FALSE)
}

load_gwas <- function(spec) {
  dt <- data.table::fread(spec$file)
  setnames(dt, old = c(
    spec$rsid_column,
    spec$effect_allele_column,
    spec$non_effect_allele_column,
    spec$effect_size_column,
    spec$standard_error_column
  ), new = c("SNP", "A1", "A2", "BETA", "SE"))
  dt[, Zscore := BETA / SE]
  dt[, .(SNP, A1, A2, BETA, SE, Zscore)]
}

# ── 2. Load GWAS files ───────────────────────────────────────────────────────
cat("[2/6] Loading GWAS files...\n")

outcome        <- load_gwas(cfg$outcome)
exposure_names <- names(cfg$exposures)
exposures      <- lapply(cfg$exposures, load_gwas)

cat(sprintf("  outcome (%s): %d SNPs\n", cfg$outcome$file, nrow(outcome)))
for (nm in exposure_names) {
  cat(sprintf("  %s (%s): %d SNPs\n", nm, cfg$exposures[[nm]]$file, nrow(exposures[[nm]])))
}

# ── 3. Harmonize alleles ─────────────────────────────────────────────────────
cat("[3/6] Harmonizing alleles against HapMap3 reference panel...\n")

data("hapmap3", package = "MRBEE", envir = environment())
ref       <- hapmap3[, c("SNP", "A1", "A2")]
gwas_list <- c(exposures, list(outcome = outcome))
gwas_list <- MRBEE::filter_align(gwas_list, ref)

n_exp   <- length(exposures)
n_snps  <- nrow(gwas_list[[1]])
cat(sprintf("  %d SNPs retained after harmonization and intersection\n", n_snps))

# ── 4. Estimate error covariance ─────────────────────────────────────────────
cat("[4/6] Estimating error covariance matrix from GWAS-insignificant SNPs...\n")

ZMatrix           <- do.call(cbind, lapply(gwas_list, `[[`, "Zscore"))
colnames(ZMatrix) <- c(exposure_names, "outcome")
Rxy               <- MRBEE::errorCov(ZMatrix)

# ── 5. IV selection via Joint.test + PLINK C+T ───────────────────────────────
cat("[5/6] Selecting instrumental variables...\n")

cat("  Running joint test across exposures...\n")
jointtest      <- MRBEE::Joint.test(
  bZ = ZMatrix[, seq_len(n_exp), drop = FALSE],
  RZ = Rxy[seq_len(n_exp), seq_len(n_exp)]
)
jointtest$SNP  <- gwas_list[[1]]$SNP

cat("  Running PLINK clump-and-threshold (p1 = 5e-8, r2 = 0.01, kb = 500)...\n")
jt_file    <- tempfile(fileext = ".txt")
clump_out  <- tempfile()
data.table::fwrite(jointtest[, c("SNP", "P")], jt_file, sep = "\t")

plink_cmd <- sprintf(
  "plink --bfile %s --clump %s --clump-field P --clump-kb 500 --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --out %s",
  cfg$ld, jt_file, clump_out
)
ret <- system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
if (ret != 0) stop("PLINK exited with an error. Check that the LD panel path is correct.",
                   call. = FALSE)

clump_file <- paste0(clump_out, ".clumped")
if (!file.exists(clump_file)) {
  stop("No SNPs passed the genome-wide significance threshold (p < 5e-8). ",
       "Check that the joint test p-values and LD panel are correct.",
       call. = FALSE)
}

ivs <- data.table::fread(clump_file)$SNP
cat(sprintf("  %d independent genome-wide significant IVs selected\n", length(ivs)))

gwas_list <- lapply(gwas_list, function(dt) dt[SNP %in% ivs])

# ── 6. MRBEE ─────────────────────────────────────────────────────────────────
cat(sprintf("[6/6] Running MRBEE.IMRP on %d IVs and %d exposure(s)...\n",
            length(ivs), n_exp))

BETA_mat <- do.call(cbind, lapply(gwas_list[seq_len(n_exp)], `[[`, "BETA"))
SE_mat   <- do.call(cbind, lapply(gwas_list[seq_len(n_exp)], `[[`, "SE"))
by       <- gwas_list[["outcome"]]$BETA
byse     <- gwas_list[["outcome"]]$SE

fit <- MRBEE::MRBEE.IMRP(
  by         = by,
  bX         = BETA_mat,
  byse       = byse,
  bXse       = SE_mat,
  Rxy        = Rxy,
  var.est    = opt[["var-est"]],
  pv.thres   = opt[["pv-thres"]],
  FDR        = !opt[["no-fdr"]]
)

results <- data.table(
  exposure  = exposure_names,
  estimate  = fit$theta,
  std_error = sqrt(diag(fit$covtheta)),
  z_stat    = fit$theta / sqrt(diag(fit$covtheta))
)
results[, p_value := 2 * pnorm(abs(z_stat), lower.tail = FALSE)]

cat("\nResults:\n")
print(results, digits = 4)

data.table::fwrite(results, opt$out)
cat("\nSaved to", opt$out, "\n")
