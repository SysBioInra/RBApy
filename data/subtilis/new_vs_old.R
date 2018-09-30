
library(tidyverse)

analyze_results <- function() {
  plot_results(read_fluxes(), read_growth_rates())
}

read_fluxes <- function() {
  fluxes <- read_reference_fluxes()
  pipeline_results <- c(
    'raw.out', 'metabolites_medium.out', 'med_met_kcat_flag.out'
  )
  for (filename in pipeline_results) {
    fluxes <- left_join(fluxes,
                        read_test_fluxes(file.path('results', filename)),
                        by='Reaction')
  }
  colnames(fluxes) <- c(
    'reaction', 'flux', 'raw', 'corrected metabolites', 'corrected kcat'
  )
  return(fluxes)
}

read_reference_fluxes <- function() {
  saturated_fluxes <- c('EdhbC', 'EmenF')
  unmatched_reactions <- c('Eatpm', 'Etrna', 'Etrna2')
  ref_fluxes <- read.table('../subtilis_ref/reactions.out', header = TRUE)
  ref_fluxes <- ref_fluxes[!(ref_fluxes$Reaction %in% unmatched_reactions), ]
  ref_fluxes[ref_fluxes$Reaction %in% saturated_fluxes, 2] <- 0
  return(ref_fluxes)
}

read_test_fluxes <- function(filename) {
  return(read.table(filename, header = TRUE))
}

read_growth_rates <- function() {
  return(c(read.table('results/growth_rates.txt')[,1]))
}

plot_results <- function(fluxes, growth_rates) {
  tidy_fluxes <- gather(fluxes, experiment, test_flux, -reaction, -flux)
  exp_levels = c('raw', 'corrected metabolites', 'corrected kcat')
  tidy_fluxes$experiment <- factor(tidy_fluxes$experiment, levels = exp_levels)
  rsquared <- sapply(exp_levels,
                     function(e) {
                       cor(tidy_fluxes$flux[tidy_fluxes$experiment == e],
                           tidy_fluxes$test_flux[tidy_fluxes$experiment == e])^2
                     })
  inset <- data.frame(experiment = exp_levels,
                      growth_rate = signif(growth_rates, 2),
                      rsquared = signif(rsquared, 2))
  result <- ggplot(tidy_fluxes) +
    geom_point(aes(flux, test_flux), color = rgb(1/255, 115/255, 178/255)) +
    geom_abline(slope=1, linetype = 2) +
    facet_wrap(vars(experiment), labeller = labeller(experiment = labels)) +
    theme_classic() +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    coord_fixed(xlim = c(-60,60)) +
    geom_text(data = inset, aes(x = -Inf, y = Inf, label = inset_label(growth_rate, rsquared)), hjust = 0, vjust = 1) +
    labs(x = 'Fluxes from hand-curated model\n(mmol/gCDW/h)',
         y = 'Fluxes from RBApy model\n(mmol/gCDW/h)')
  pdf('results/pipeline_vs_ref.pdf', 9, 4); print(result); dev.off()
}

inset_label <- function(growth_rate, rsquared) {
  return(lapply(seq_along(growth_rate), function(i) paste0(
    "growth rate = ", growth_rate[i], " 1/h\n", "R^2 = ", rsquared[i]
  )))
}

analyze_results()
