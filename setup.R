# Explicit setup entry point; never runs automatically when sourcing model code.
# Non-interactive Rscript cannot answer a prompt to create a personal library.
if (file.access(.libPaths()[1], 2) != 0) {
  user_library <- path.expand(Sys.getenv("R_LIBS_USER"))
  if (!nzchar(user_library)) stop("Set R_LIBS_USER to a writable R library directory")
  dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_library, .libPaths()))
}
message("R library: ", .libPaths()[1])
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("posterior", quietly = TRUE)) install.packages("posterior")
if (!requireNamespace("cmdstanr", quietly = TRUE))
  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::check_cmdstan_toolchain()
if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  cmdstanr::install_cmdstan(cores = 2)

message("Setup complete. CmdStan: ", cmdstanr::cmdstan_path())
message("Next: Rscript tests.R, then Rscript examples/one_fit.R")
