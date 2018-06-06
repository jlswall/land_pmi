library("tidyverse")


## When I organized the data, any taxa which didn't get to at least
## at certain frequency (I made it 1%) for at least one day and one
## individual during *the whole time period* got folded into the
## "Rare" taxa group.  Now, I'm only using data for the first couple
## of weeks, so I need to check that all of these taxa have at least
## one day during the first couple of weeks where the frequency
## exceeds the cutoff for at least one individual.  Otherwise, those
## taxa should be considered "rare" for an analysis of just the first
## 2 weeks.
freqCutoff <- 0.01



## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "orders"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0("../../../", taxalevel, "_massaged.csv"))

## We're only looking at the first 15 days (approx. 2 weeks).
earlyT <- taxaT %>% filter(days <= 15)

rm(taxaT)
## ##################################################



## ##################################################
## Any taxa which don't make up at least a certain percentage (given
## by freqCutoff) on at least 1 day and 1 cadaver (during the first 2
## weeks) will get absorbed into the "rare" category.

## Get list of taxa for which the percentage never exceeds the
## frequency cutoff, for any day or individual.
freqTaxaT <- earlyT %>%
  group_by(taxa) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  filter(maxFracBySubjDay >= freqCutoff) %>%
  arrange(desc(maxFracBySubjDay)) %>%
  select(taxa)


## Incorporate the taxa that occur less than the frequency cutoff into
## the "rare" group.
commontaxaT <- earlyT
commontaxaT[!(commontaxaT$taxa %in% freqTaxaT$taxa), "taxa"] <- "Rare"
commontaxaT <- commontaxaT %>%
  group_by(days, degdays, subj, taxa) %>%
  summarize(counts=sum(counts), fracBySubjDay=sum(fracBySubjDay)) %>%
  ungroup()


## Check that the fractions add up to 1, appropriately.
unique(
    commontaxaT %>%
    group_by(days, subj) %>%
    summarize(sumFracBySubjDay = sum(fracBySubjDay)) %>%
    ungroup() %>%
    select(sumFracBySubjDay)
)


## Remove the list of taxa names that satisfied the frequency cutoff.
rm(freqTaxaT, earlyT)
## ##################################################



## ##################################################
## Write out the data to use in analyses for the first 2 weeks.

write.csv(commontaxaT, file="orders_only_first_two_weeks.csv", row.names=FALSE)
## ##################################################
