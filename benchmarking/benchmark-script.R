library(coiaf)
library(McCOILR)

# Submit job to slurm
benchmark <- function(n_loci, n_samples, min_iterations) {
  set.seed(500)

  missingness <- 0.05
  coi <- as.integer(runif(n_samples, 1, 10))
  plmaf <- runif(n_loci, 0, 0.5)

  # create genotype matrix for this using sim biallelic
  f <- lapply(seq_len(n_samples), function(x) {
    coiaf::sim_biallelic(coi[x], plmaf, epsilon = 0.05)
  })

  # now let's create out genotype matrix for RMCL
  gts <- lapply(f, function(x) {
    gtmat <- x$data$wsmaf
    gtmat[gtmat >= 0.95] <- 1
    gtmat[gtmat <= 0.05] <- 1
    gtmat[gtmat < 0.95 & gtmat > 0.05] <- 0.5
    gtmat[which(as.logical(rbinom(100, 1, 0.05)))] <- -1
    return(gtmat)
  })

  # group together to get out gtmat
  gtmat <- do.call(rbind, gts)
  rownames(gtmat) <- letters[seq_len(nrow(gtmat))]
  colnames(gtmat) <- paste0("pos_", seq_len(ncol(gtmat)))
  gtmat <- as.data.frame(gtmat)


  bench::mark(
    rmcl = McCOIL_categorical(
      data = gtmat,
      maxCOI = 25,
      totalrun = 100000,
      burnin = 1000,
      err_method = 3,
      path = tempdir(),
    ),
    coiaf_disc = lapply(f, coiaf::compute_coi, "sim"),
    coiaf_cont = lapply(f, coiaf::optimize_coi, "sim"),
    min_iterations = min_iterations,
    check = FALSE
  )
}

results <- benchmark(n_loci = 1000, n_samples = 100, min_iterations = 10)
