## compare_branches.R
## Sources the R/ files from each branch into isolated environments and
## checks that every exported function returns numerically identical results.
##
## Usage:
##   Rscript compare_branches.R [options]
##
## Prerequisites: MASS (already in renv.lock), optparse (install separately)

library(MASS)
library(optparse)

opt_list <- list(
  make_option(c("-a", "--branch-a"), type = "character", default = "HEAD",
              metavar = "BRANCH",
              help    = "First git branch or ref  [default: %default]"),
  make_option(c("-b", "--branch-b"), type = "character", default = "main",
              metavar = "BRANCH",
              help    = "Second git branch or ref [default: %default]"),
  make_option(c("--label-a"),        type = "character", default = NULL,
              metavar = "LABEL",
              help    = "Display name for branch A [default: -a value]"),
  make_option(c("--label-b"),        type = "character", default = NULL,
              metavar = "LABEL",
              help    = "Display name for branch B [default: -b value]"),
  make_option(c("-t", "--tol"),      type = "double",    default = 1e-8,
              metavar = "TOL",
              help    = "Numerical tolerance for all.equal() [default: %default]")
)

opts <- parse_args(OptionParser(option_list = opt_list))

repo_root <- trimws(system("git rev-parse --show-toplevel", intern = TRUE))
worktrees <- character(0)

on.exit({
  for (wt in worktrees) {
    system2("git", c("-C", repo_root, "worktree", "remove", "--force", wt),
            stdout = FALSE, stderr = FALSE)
  }
}, add = TRUE)

# Creates a worktree for a git ref and returns its path.
resolve_branch <- function(arg) {
  # Use the SHA to avoid "already checked out" errors.
  sha <- tryCatch(
    trimws(system2("git", c("-C", repo_root, "rev-parse", "--verify", arg),
                   stdout = TRUE, stderr = TRUE)),
    warning = function(w) character(0)
  )
  if (length(sha) == 0 || !grepl("^[0-9a-f]{40}$", sha)) {
    stop("'", arg, "' is not a known git ref.")
  }

  wt_path <- file.path(tempdir(), paste0("mrbee_cmp_", gsub("[^A-Za-z0-9_.-]", "_", arg)))
  ret <- system2("git", c("-C", repo_root, "worktree", "add", "--detach", wt_path, sha),
                 stdout = TRUE, stderr = TRUE)
  if (!dir.exists(wt_path)) {
    stop("git worktree add failed for '", arg, "':\n", paste(ret, collapse = "\n"))
  }
  worktrees <<- c(worktrees, wt_path)
  wt_path
}

path_a  <- resolve_branch(opts[["branch-a"]])
path_b  <- resolve_branch(opts[["branch-b"]])
label_a <- if (!is.null(opts[["label-a"]])) opts[["label-a"]] else opts[["branch-a"]]
label_b <- if (!is.null(opts[["label-b"]])) opts[["label-b"]] else opts[["branch-b"]]
tol     <- opts[["tol"]]

source_pkg <- function(path) {
  env <- new.env(parent = baseenv())
  env$library <- function(...) invisible(NULL)   # suppress re-loading MASS
  for (f in list.files(file.path(path, "R"), pattern = "\\.R$", full.names = TRUE)) {
    sys.source(f, envir = env, chdir = FALSE)
  }
  env
}

message("Sourcing ", label_a, " (", path_a, ") …")
env_admin <- source_pkg(path_a)

message("Sourcing ", label_b, " (", path_b, ") …")
env_main <- source_pkg(path_b)

# ── shared synthetic data ────────────────────────────────────────────────────

set.seed(1234)
n <- 200
p <- 3

bX         <- matrix(rnorm(n * p, sd = 0.5), n, p)
bXse       <- matrix(0.1, n, p)
theta_true <- c(0.3, -0.2, 0.1)
by         <- drop(bX %*% theta_true) + rnorm(n, sd = 0.1)
byse       <- rep(0.1, n)
Rxy        <- diag(p + 1)

bx_uv   <- bX[, 1]
bxse_uv <- bXse[, 1]
Rxy_uv  <- diag(2)

# ── comparison helper ────────────────────────────────────────────────────────

compare <- function(label, a, b) {
  chk <- all.equal(a, b, tolerance = tol, check.attributes = FALSE)
  ok  <- isTRUE(chk)
  if (ok) message("  PASS  ", label) else { message("  FAIL  ", label); message(chk) }
  ok
}

results <- list()

# ── MRBEE.IMRP ───────────────────────────────────────────────────────────────

message("\n── MRBEE.IMRP ──────────────────────────────────────────────────────────────")
args_mv <- list(by = by, bX = bX, byse = byse, bXse = bXse, Rxy = Rxy)

res_a <- do.call(env_admin$MRBEE.IMRP, args_mv)
res_m <- do.call(env_main$MRBEE.IMRP,  args_mv)

results[["MRBEE.IMRP / theta"]]    <- compare("theta",    res_a$theta,    res_m$theta)
results[["MRBEE.IMRP / covtheta"]] <- compare("covtheta", res_a$covtheta, res_m$covtheta)
results[["MRBEE.IMRP / delta"]]    <- compare("delta",    res_a$delta,    res_m$delta)

# ── MRBEE.IMRP.Egger ─────────────────────────────────────────────────────────

message("\n── MRBEE.IMRP.Egger ────────────────────────────────────────────────────────")

res_ae <- do.call(env_admin$MRBEE.IMRP.Egger, args_mv)
res_me <- do.call(env_main$MRBEE.IMRP.Egger,  args_mv)

results[["MRBEE.IMRP.Egger / theta"]]    <- compare("theta",    res_ae$theta,    res_me$theta)
results[["MRBEE.IMRP.Egger / covtheta"]] <- compare("covtheta", res_ae$covtheta, res_me$covtheta)
results[["MRBEE.IMRP.Egger / delta"]]    <- compare("delta",    res_ae$delta,    res_me$delta)

# ── MRBEE.IMRP.UV ────────────────────────────────────────────────────────────

message("\n── MRBEE.IMRP.UV ───────────────────────────────────────────────────────────")
args_uv <- list(by = by, bx = bx_uv, byse = byse, bxse = bxse_uv, Rxy = Rxy_uv)

res_au <- do.call(env_admin$MRBEE.IMRP.UV, args_uv)
res_mu <- do.call(env_main$MRBEE.IMRP.UV,  args_uv)

results[["MRBEE.IMRP.UV / theta"]]    <- compare("theta",    res_au$theta,    res_mu$theta)
results[["MRBEE.IMRP.UV / vartheta"]] <- compare("vartheta", res_au$vartheta, res_mu$vartheta)
results[["MRBEE.IMRP.UV / delta"]]    <- compare("delta",    res_au$delta,    res_mu$delta)

# ── summary ──────────────────────────────────────────────────────────────────

message("\n── Summary ─────────────────────────────────────────────────────────────────")
n_pass <- sum(unlist(results))
n_fail <- sum(!unlist(results))
message(sprintf("  %d / %d checks passed", n_pass, n_pass + n_fail))

if (n_fail > 0) {
  message("\nFailed checks:")
  for (nm in names(results)[!unlist(results)]) message("  - ", nm)

  message("\n", label_a, " theta: ", paste(round(res_a$theta, 6), collapse = ", "))
  message(label_b, " theta: ", paste(round(res_m$theta, 6), collapse = ", "))
  quit(status = 1)
} else {
  message("\nBoth branches produce identical results.")
}
