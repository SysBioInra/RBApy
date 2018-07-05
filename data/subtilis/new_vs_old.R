
saturated_fluxes <- c('EdhbC', 'EmenF')
unmatched_reactions <- c('Eatpm', 'Etrna', 'Etrna2', 'R_maintenance_atp')

ref_fluxes <- read.table('../subtilis_ref/reactions.out', header = TRUE)
auto_fluxes <- read.table('reactions.out', header = TRUE)

# remove reactions that have no match
ref_fluxes <- ref_fluxes[!(ref_fluxes$Reaction %in% unmatched_reactions), ]
auto_fluxes <- auto_fluxes[!(auto_fluxes$Reaction %in% unmatched_reactions), ]

# set arbitrarily saturated fluxes to 0 in reference model
ref_fluxes[ref_fluxes$Reaction %in% saturated_fluxes, 2] <- 0

# match reactions
reaction_match <- match(ref_fluxes$Reaction, auto_fluxes$Reaction)
auto_fluxes <- auto_fluxes[reaction_match, ]

# linear model
model <- lm(auto_fluxes$Flux ~ ref_fluxes$Flux)

# plot results
pdf('pipeline_vs_ref.pdf')
plot(ref_fluxes$Flux, auto_fluxes$Flux)
abline(model)
dev.off()
