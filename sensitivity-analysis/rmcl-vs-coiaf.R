library(coiaf)

# Define functions -------------------------------------------------------------
sa <- function(coverage, loci, overdispersion, relatedness, seq_error) {
  param_grid <- tidyr::expand_grid(
    coverage = coverage,
    loci = loci,
    overdispersion = overdispersion,
    relatedness = relatedness,
    seq_error = seq_error
  )

  purrr::pmap(param_grid, run_estimations)
}

run_estimations <- function(coverage,
                            loci,
                            overdispersion,
                            relatedness,
                            seq_error) {
  rlang::check_required(coverage)
  rlang::check_required(loci)
  rlang::check_required(overdispersion)
  rlang::check_required(relatedness)
  rlang::check_required(seq_error)

  withr::local_seed(2022)

  # Setup inputs
  n_cois <- 200
  coi <- as.integer(rep(1:20, 10))
  coi <- rlang::set_names(coi, ~ paste0("sample_", seq_along(.)))
  plmaf <- runif(loci, 0, 0.5)

  # Simulate data
  sim_data <- purrr::map(
    coi,
    coiaf::sim_biallelic,
    plmaf = plmaf,
    coverage = coverage,
    overdispersion = overdispersion,
    relatedness = relatedness,
    epsilon = seq_error
  )

  # Compute genotypes
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

  # Estimated plmaf
  coiaf_plmaf <- sim_data %>%
    purrr::map(~ .x$data$wsmaf) %>%
    do.call(rbind, .)

  # Estimate seq error on each sample using
  estimated_seq_error <- purrr::map_dbl(sim_data, function(x) {
    coiaf:::estimate_seq_error(
      x$data$wsmaf,
      colMeans(coiaf_plmaf),
      bin_size = 0.005 * loci
    )
  }) %>%
    median()

  # Run coiaf
  coiaf_disc <- purrr::map(
    sim_data,
    coiaf::compute_coi,
    data_type = "sim",
    seq_error = estimated_seq_error
  )
  coiaf_cont <- purrr::map(
    sim_data,
    coiaf::optimize_coi,
    data_type = "sim",
    seq_error = estimated_seq_error
  )

  # Run RMCL
  rmcl <- McCOILR::McCOIL_categorical(
    data = gtmat,
    maxCOI = 25,
    totalrun = 5000,
    burnin = 2000,
    M0 = 5,
    thin = 0.1,
    err_method = 3,
    path = tempdir()
  )

  # Get coi estimates
  coiaf_disc_res <- coiaf_disc %>%
    purrr::imap(~ .x[[1]]) %>%
    purrr::flatten_dbl()
  coiaf_cont_res <- purrr::flatten_dbl(coiaf_cont)
  rmcl_res <- rmcl[seq_len(n_cois),]

  # Combine results into tibble
  coi_res <- rmcl_res %>%
    dplyr::select(-c(file, CorP)) %>%
    tibble::add_column(
      coiaf_disc = coiaf_disc_res,
      coiaf_cont = coiaf_cont_res,
      truth = coi
    )

  # Adjust estimated plmafs
  coiaf_plamf_adj <- colMeans(coiaf_plmaf) -
    (((0.5 - colMeans(coiaf_plmaf)) * 2) * seq_error)

  plmaf_res <- tibble::tibble(
    loci = rmcl$name[rmcl$CorP == "P"],
    RMCL = rmcl$median[rmcl$CorP == "P"]
  ) %>%
    dplyr::mutate(
      loci_n = stringr::str_extract(loci, "(?<=_)\\d+"),
      loci_n = as.double(loci_n)
    ) %>%
    dplyr::left_join(data.frame("loci_n" = seq_len(loci), "truth" = plmaf)) %>%
    dplyr::left_join(data.frame("loci_n" = seq_len(loci), "coiaf" = coiaf_plamf_adj))

  list(coi_res, plmaf_res)
}

plot_cois <- function(data, title = "") {
  data %>%
    dplyr::rename(sample = name, RMCL = median) %>%
    tidyr::pivot_longer(cols = c(coiaf_disc, coiaf_cont, RMCL)) %>%
    dplyr::group_by(truth, name) %>%
    dplyr::summarise(
      m = median(value),
      q1 = quantile(value, 0.025),
      q2 = quantile(value, 0.975)
    ) %>%
    ggplot2::ggplot(ggplot2::aes(truth, m, ymin = q1, ymax = q2, color = name)) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.1)) +
    ggplot2::scale_colour_discrete(
      breaks = c("coiaf_disc", "coiaf_cont", "RMCL"),
      labels = c(
        expression("Discrete"~italic("coiaf")),
        expression("Continuous"~italic("coiaf")),
        expression(italic("THE REAL McCOIL"))
      )
      ) +
    ggplot2::labs(x = "True COI", y = "Estimated COI", title = title) +
    theme_coiaf() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
}

# Run the sensitivity analysis -------------------------------------------------
results <- sa(
  coverage = 200,
  # coverage = c(50, 100, 500, 1000),
  loci = 1000,
  # loci = c(100, 1000, 10000),
  overdispersion = 0,
  # overdispersion = c(0.01, 0.05, 0.1, 0.15),
  relatedness = 0,
  # relatedness = c(0.05, 0.25, 0.5),
  seq_error = 0
  # seq_error = c(0.005, 0.01, 0.015)
)

saveRDS(
  results,
  here::here(
    "sensitivity-analysis",
    "sa-results",
    glue::glue("seqerr.rds")
  )
)

# Make figures -----------------------------------------------------------------
# Generate figures for SA
seqerr %>%
  purrr::imap(~ .x[[1]]) %>%
  purrr::map2(paste("Sequence Error =", c(0.005, 0.01, 0.015)), plot_cois) %>%
  patchwork::wrap_plots(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = 'bold')
  )

ggplot2::ggsave(
  filename = here::here("sensitivity-analysis", "sa-results", "seqerr.png"),
  device = "png",
  width = 2500,
  height = 1500,
  units = "px",
  dpi = "print"
)

# Plot PLMAF coiaf vs RMCL
seqerr %>%
  purrr::imap(~ .x[[2]]) %>%
  purrr::map2(paste("Sequence Error =", c(0.005, 0.01, 0.015)), function(x, title) {
    x %>%
      tidyr::pivot_longer(
        cols = c("RMCL", "coiaf"),
        names_to = "package",
        values_to = "estimate"
      ) %>%
      dplyr::mutate(
        package = dplyr::recode(
          package,
          coiaf = "coiaf",
          RMCL = "THE REAL McCOIL"
        )
      ) %>%
      ggplot(aes(truth, estimate, color = package, group = package)) +
      geom_point(alpha = 0.25) +
      geom_smooth(method = "lm", se = FALSE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      facet_wrap(~package, ncol = 2) +
      guides(color = "none") +
      labs(x = "True PLMAF", y = "Estimated PLMAF", title = title) +
      theme_coiaf() +
      theme(legend.title = element_blank())
  }) %>%
  patchwork::wrap_plots(nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = 'bold')
  )

ggplot2::ggsave(
  filename = here::here(
    "sensitivity-analysis",
    "sa-results",
    "seqerr-plmaf.png"
  ),
  device = "png",
  width = 2500,
  height = 1500,
  units = "px",
  dpi = "print"
)
