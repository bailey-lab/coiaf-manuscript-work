library(coiaf)
library(McCOILR)

# Declare file location
here::i_am("benchmarking/benchmark-script.R")

# Initialize variables
# Note that these variables can be controlled from the global environment if
# the user chooses to run an RStudio job
# n_loci = 100
# n_samples = 100

# Print info on variables
cli::cli_inform("Running {n_loci} loci and {n_samples} samples.")

# Submit job to slurm
benchmark <- function(n_loci, n_samples, min_iterations) {
  withr::local_output_sink(new = nullfile())
  withr::local_seed(500)

  coi <- as.integer(runif(n_samples, 1, 10))
  coi <- rlang::set_names(coi, ~ paste0("sample_", seq_along(.)))
  plmaf <- runif(n_loci, 0, 0.5)

  # create genotype matrix for this using sim biallelic
  sim_data <- purrr::map(coi, ~ coiaf::sim_biallelic(.x, plmaf))

  # Compute genotypes
  seq_error <- 0.01
  genotypes <- purrr::map(sim_data, function(x) {
    gtmat <- x$data$wsmaf
    gtmat[gtmat >= (1 - seq_error)] <- 1
    gtmat[gtmat <= seq_error] <- 0
    gtmat[gtmat > seq_error & gtmat < (1 - seq_error)] <- 0.5
    gtmat
  })

  # Get genotype matrix to feed into RMCL
  gtmat <- do.call(rbind, genotypes)
  colnames(gtmat) <- paste0("pos_", seq_len(ncol(gtmat)))
  gtmat <- as.data.frame(gtmat)

  bench::mark(
    rmcl = McCOIL_categorical(
      data = gtmat,
      maxCOI = 25,
      totalrun = 5000,
      burnin = 2000,
      M0 = 5,
      err_method = 3,
      path = tempdir()
    ),
    coiaf_disc = lapply(sim_data, coiaf::compute_coi, data_type = "sim"),
    coiaf_cont = lapply(sim_data, coiaf::optimize_coi, data_type = "sim"),
    coiaf_boot = lapply(sim_data, coiaf::bootstrap_ci, parallel = FALSE),
    coiaf_boot_par = lapply(sim_data, coiaf::bootstrap_ci, parallel = TRUE),
    min_iterations = min_iterations,
    check = FALSE,
    memory = FALSE
  )
}

results <- benchmark(
  n_loci = n_loci,
  n_samples = n_samples,
  min_iterations = 10
)

# Save data
saveRDS(
  results,
  here::here(
    "benchmarking",
    "results",
    glue::glue("l{ n_loci }_s{ n_samples }.rds")
  )
)
