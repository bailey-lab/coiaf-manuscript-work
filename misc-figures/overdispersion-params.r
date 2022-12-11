# Plot for impact of changing alpha parameter ----------------------------------
LaplacesDemon::ddirichlet(c(.1, .3, .6), rep(1, 3))
LaplacesDemon::rdirichlet(10, rep(1, 3))

# Set params
cois <- c(2, 3, 5, 10)
alphas <- c(0.01, 0.025, 0.05, 0.1, 0.5, 1, 5, 10, 20, 50, 100, 1000)

pars <- expand.grid("coi" = cois, "alpha" = alphas)
pars$mean <- 0
pars$sd <- 0

# Generate data
library(coiaf)
for (i in seq_along(pars$coi)) {
  m <- replicate(1000, max(coiaf:::rdirichlet(rep(pars$alpha[i], pars$coi[i]))))
  pars$mean[i] <- mean(m)
  pars$sd[i] <- sd(m)
}

# Plot
library(tidyverse)
ggplot(pars, aes(alpha, mean, color = as.factor(coi))) +
  geom_hline(
    data = data.frame(coi = cois, y = 1 / cois),
    aes(yintercept = y, color = as.factor(coi)),
    linetype = "dashed"
  ) +
  geom_line() +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks(sides = "b") +
  ylab("Proportion of total parasite density contributed \nby highest density parasite strain") +
  # xlab("Alpha (Parameter controlling overdispersion \nof parasite densities)") +
  xlab("Alpha") +
  # ggtitle("Dashed lines show maximum parasite density if \nall strains had the same parasite density, i.e. 1/COI") +
  scale_color_discrete(name = "COI") +
  theme_coiaf() +
  theme(plot.title = element_text(size = 10))

# Save
ggsave(
  filename = here::here("misc-figures", "alpha-param.png"),
  device = "png",
  width = 2000,
  height = 1750,
  units = "px",
  dpi = "print"
)

# Plot for impact of changing overdispersion parameter -------------------------
# Set params and get data
overd <- c(0.5, 1, 5, 10, 100)
counts <- 0:100
res <- list()
for (i in seq_along(overd)) {
  res[[i]] <- data.frame(
    "prob" = emdbook::dbetabinom(
      counts,
      prob = 0.5,
      theta = overd[i],
      size = max(counts)
    ),
    "counts" = counts,
    "overdispersion" = overd[i]
  )
}

# Plot
do.call(rbind, res) %>%
  ggplot(aes(counts, prob, color = as.factor(overdispersion))) +
  geom_vline(xintercept = 50, linetype = "dashed") +
  geom_line() +
  ylab("Probability Density") +
  xlab("Reference Reads") +
  scale_color_discrete(name = "Overdispersion\nin Sequencing") +
  theme_coiaf() +
  ggtitle("Sequence Coverage = 100, Expected Reference Reads = 50") +
  theme(plot.title = element_text(size = 10))

# Save
ggsave(
  filename = here::here("misc-figures", "over-param.png"),
  device = "png",
  width = 2000,
  height = 1750,
  units = "px",
  dpi = "print"
)
