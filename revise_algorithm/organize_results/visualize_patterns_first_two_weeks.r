library("tidyverse")

## ##################################################
## Combined taxa.

## Read in the combined data in wide format.
widetaxaT <- read_csv("../all_together/combined_taxa.csv")

## For ease of making graphs, put in long format.  Also, subset to the
## first 15 days.
earlyT <- gather(widetaxaT, taxa, fracBySubjDay, -days, -degdays, -subj) %>% filter(days<=15)

rm(widetaxaT)
## ##################################################


## ##################################################
## Read list of taxa and groups they belong to (orders, phyla, etc.)

taxaGrpsT <- read_csv("../all_together/list_taxas_groups.csv")
## ##################################################



## ##################################################
## Make plots.
## For the first 15 days, using orders gave the best percentage of
## explained variation.


## ###################
## For original units, the most important taxa seem to be:
impCombTaxa <- c("JG30_KF_CM45", "Rhizobiales", "Actinomycetales", "Pseudomonadales", "Bacteroidales")
## Subset to just these taxa.
graphT <- earlyT %>% filter(taxa %in% impCombTaxa)
## Add in the group associated with the taxa.
graphT <- graphT %>%
  left_join(taxaGrpsT)

## Make graph.
ggplot(graphT, aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=taxa)) +
##  facet_grid(~group) + 
  xlab("Degree days") +
  ylab("Fraction by degree day and subject") +
  ggtitle("Top 5 taxa using combined taxa (first 15 days, orig. units)")
## Save to a PDF file.
ggsave("orders_orig_units_first_two_weeks_top_5_taxa.pdf", width=7, height=3.5, units="in")
rm(graphT, impCombTaxa)
## ###################


## ###################
## For square rooot units, the most important taxa seem to be:
impCombTaxa <- c("JG30_KF_CM45", "Bacteroidales", "Burkholderiales", "Rhizobiales", "Actinomycetales")
## Subset to just these taxa.
graphT <- earlyT %>% filter(taxa %in% impCombTaxa)
## Add in the group associated with the taxa.
graphT <- graphT %>%
  left_join(taxaGrpsT)

## Make graph.
ggplot(graphT, aes(sqrt(degdays), fracBySubjDay)) +
  geom_point(aes(color=taxa)) +
##  facet_wrap(~group, scales="free_y") + 
  xlab("Square root of degree days") +
  ylab("Fraction by degree day and subject") +
  ggtitle("Top 5 taxa using combined taxa (first 15 days, sqrt. units)")
## Save to a PDF file.
ggsave("orders_sqrt_units_first_two_weeks_top_5_taxa.pdf", width=7, height=3.5, units="in")
rm(graphT, impCombTaxa)
## ##################################################
