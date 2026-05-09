lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
stopifnot(nchar(lib) > 1)

install.packages(
  c("FDRestimation"),
  repos = "https://cloud.r-project.org",
  lib = lib
)

install.packages(
  "../",
  repos = NULL,
  type = "source",
  lib = lib
)

cat("Setup complete.\n")
