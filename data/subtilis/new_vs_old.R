
library(tidyverse)

analyze_results <- function() {
  plot_results(read_fluxes(), read_growth_rates())
}

read_fluxes <- function() {
  fluxes <- read_reference_fluxes()
  pipeline_results <- c(
    'raw.out', 'metabolites_medium.out', 'tmp.out'
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
                      label = sprintf('growth rate = %.2f\nR^2 = %.2f',
                                      growth_rates, rsquared))
  labels <- c('raw' = 'First run',
              'corrected metabolites' = 'After matching of metabolite IDs',
              'corrected kcat' = 'After adjustment of catalytic constants')
  result <- ggplot(tidy_fluxes) +
    geom_point(aes(flux, test_flux), color = 'blue') +
    geom_abline(slope=1) +
    facet_wrap(vars(experiment), labeller = labeller(experiment = labels)) +
    geom_text(data = inset, aes(x = -38, y = 60, label = label), hjust = 0, vjust = 1) +
    labs(x = 'Fluxes from hand-curated model\n(mmol/gCDW/h)',
         y = 'Fluxes from RBApy model\n(mmol/gCDW/h)')
  pdf('results/pipeline_vs_ref.pdf', 9, 4); print(result); dev.off()
}

analyze_results()
