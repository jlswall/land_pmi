library("tidyverse")
library("randomForest")
library("figdim")


## ##################################################
## Read in different levels of taxa.

## #####################
## Read in phyla-level taxa.
taxalevel <- "phyla"
## Read in cleaned-up phyla, orders, or families taxa.
tmpT <- read_csv(paste0("../", taxalevel, "_massaged.csv"))
## Move to wide format.
phylaT <- tmpT %>%
  select(days, degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay)
## Also, save the list of unique phyla taxa as a data frame.
taxasGroupT <- tmpT %>% distinct(taxa)
taxasGroupT$group <- "phylum"
## #####################

## #####################
## Read in order-level taxa.
taxalevel <- "orders"
## Read in cleaned-up phyla, orders, or families taxa.
tmpT <- read_csv(paste0("../", taxalevel, "_massaged.csv"))
## Move to wide format.
ordersT <- tmpT %>%
  select(days, degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay)
## Save the list of unique order taxa as a tibble, combine it with
## those of the phyla.
tmpGroupT <- tmpT %>% distinct(taxa)
tmpGroupT$group <- "order"
taxasGroupT <- bind_rows(taxasGroupT, tmpGroupT)
rm(tmpGroupT)
## #####################

## #####################
## Read in order-level taxa.
taxalevel <- "families"
## Read in cleaned-up phyla, orders, or families taxa.
tmpT <- read_csv(paste0("../", taxalevel, "_massaged.csv"))
## Move to wide format.
familiesT <- tmpT %>%
  select(days, degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay)
## Save the list of unique family taxa as a tibble, combine it with
## those of the phyla and orders.
tmpGroupT <- tmpT %>% distinct(taxa)
tmpGroupT$group <- "family"
taxasGroupT <- bind_rows(taxasGroupT, tmpGroupT)
rm(tmpGroupT)
## #####################

rm(taxalevel, tmpT)
## ##################################################


## ##################################################
## Write out the list of taxas with their group.

write.csv(taxasGroupT %>% filter(taxa!="Rare"), file="list_taxas_groups.csv", row.names=FALSE)
## ##################################################



## ##################################################
## Combine all the taxa.

## Merging on days, degdays, and subj.  Remove the "Rare" category
## before doing the merge.
allT <- phylaT %>%
  select(-Rare) %>%
  inner_join(ordersT %>% select(-Rare), by=c("days", "degdays", "subj")) %>%
  inner_join(familiesT %>% select(-Rare), by=c("days", "degdays", "subj"))
## ##################################################



write.csv(allT, file="combined_taxa.csv", row.names=FALSE)

