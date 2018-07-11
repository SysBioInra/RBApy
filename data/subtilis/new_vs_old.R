
library(tidyverse)

analyze_results <- function() {
  plot_results(read_fluxes())
}

read_fluxes <- function() {
  fluxes <- read_reference_fluxes()
  pipeline_results <- c(
    'raw.out', 'metabolites_medium.out', 'med_met_kcat_flag.tsv'
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

plot_results <- function(fluxes) {
  tidy_fluxes <- gather(fluxes, experiment, test_flux, -reaction, -flux)
  tidy_fluxes$experiment <- factor(
    tidy_fluxes$experiment,
    levels=c('raw', 'corrected metabolites', 'corrected kcat'))
  pdf('results/pipeline_vs_ref.pdf')
  print(ggplot(tidy_fluxes) + geom_point(aes(flux, test_flux))
        + geom_abline(slope=1) +
        facet_wrap(vars(experiment))
        )
  dev.off()
}

analyze_results()
