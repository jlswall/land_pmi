library("tidyverse")

## ##################################################
## Combined taxa.

## Read in the combined data in wide format.
widetaxaT <- read_csv("../all_together/combined_taxa.csv")

## For ease of making graphs, put in long format.
taxaT <- gather(widetaxaT, taxa, fracBySubjDay, -days, -degdays, -subj)

rm(widetaxaT)
## ##################################################


## ##################################################
## Read list of taxa and groups they belong to (orders, phyla, etc.)

taxaGrpsT <- read_csv("../all_together/list_taxas_groups.csv")
## ##################################################



## ##################################################
## Make plots.

## ###################
## For original units, the most important taxa seem to be:
impCombTaxa <- c("Listeriaceae", "Peptostreptococcaceae", "Dietziaceae", "Chloroflexi", "Microbacteriaceae")
## Subset to just these taxa.
graphT <- taxaT %>% filter(taxa %in% impCombTaxa)
## Add in the group associated with the taxa.
graphT <- graphT %>%
  left_join(taxaGrpsT)

## Make graph.
ggplot(graphT, aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=taxa)) +
  facet_grid(~group) + 
  xlab("Degree days") +
  ylab("Fraction by degree day and subject") +
  ggtitle("Top 5 taxa using combined taxa (all time steps, orig. units)")
## Save to a PDF file.
ggsave("combined_orig_units_all_data_top_5_taxa.pdf", width=7, height=3.5, units="in")
rm(graphT, impCombTaxa)
## ###################


## ###################
## For square rooot units, the most important taxa seem to be:
impCombTaxa <- c("Listeriaceae", "Peptostreptococcaceae", "Dietziaceae", "Firmicutes", "Turicibacterales")
## Subset to just these taxa.
graphT <- taxaT %>% filter(taxa %in% impCombTaxa)
## Add in the group associated with the taxa.
graphT <- graphT %>%
  left_join(taxaGrpsT)

## Make graph.
ggplot(graphT, aes(sqrt(degdays), fracBySubjDay)) +
  geom_point(aes(color=taxa)) +
  facet_wrap(~group, scales="free_y") + 
  xlab("Square root of degree days") +
  ylab("Fraction by degree day and subject") +
  ggtitle("Top 5 taxa using combined taxa (all time steps, sqrt. units)")
## Save to a PDF file.
ggsave("combined_sqrt_units_all_data_top_5_taxa.pdf", width=7, height=3.5, units="in")
rm(graphT, impCombTaxa)




impCombTaxa <- c("Listeriaceae", "Peptostreptococcaceae", "Turicibacterales", "JG30_KF_CM45", "Firmicutes", "Chloroflexi")

## ##################################################
