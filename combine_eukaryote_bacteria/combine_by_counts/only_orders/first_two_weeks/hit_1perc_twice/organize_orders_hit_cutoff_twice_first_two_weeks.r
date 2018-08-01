library("tidyverse")


## I'm only using data for the first couple of weeks, so I need to
## check that, for all of these taxa, each exceeds the frequency
## cutoff for at least two separate subjects (can be same or different
## days).  Otherwise, those taxa should be considered "rare" for this
## analysis.
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
## by freqCutoff) on at least 1 day and 2 cadavers (during the first 2
## weeks) will get absorbed into the "rare" category.

## Get list of taxa for which the percentage exceeds the frequency
## cutoff, for at least 2 individuals.
freqTaxaT <- earlyT %>%
  group_by(taxa, subj) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  group_by(taxa) %>%
  summarize(ctSubjExceed = sum(maxFracBySubjDay > freqCutoff)) %>%
  filter(ctSubjExceed > 1) %>% ## arrange(desc(ctSubjExceed)) %>%
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

write.csv(commontaxaT, file="orders_hit_cutoff_twice_first_two_weeks.csv", row.names=FALSE)
## ##################################################
