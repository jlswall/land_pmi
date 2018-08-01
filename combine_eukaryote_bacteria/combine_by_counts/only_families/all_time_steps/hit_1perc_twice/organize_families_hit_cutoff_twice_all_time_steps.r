library("tidyverse")


## I need to check that all the taxa have at least one day for at
## least TWO individuals on which the frequency exceeds the cutoff.
## Otherwise, those taxa should be grouped as "rare".
freqCutoff <- 0.01



## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0("../../../", taxalevel, "_massaged.csv"))
## ##################################################



## ##################################################
## Any taxa which don't make up at least a certain percentage (given
## by freqCutoff) on at least 1 day for at least TWO cadavers will get
## absorbed into the "rare" category.

## Get list of taxa for which the percentage exceeds the
## frequency cutoff for at least 2 individuals.
freqTaxaT <- taxaT %>%
  group_by(taxa, subj) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  group_by(taxa) %>%
  summarize(ctSubjExceed = sum(maxFracBySubjDay > freqCutoff)) %>%
  filter(ctSubjExceed > 1) %>%
  select(taxa)


## Incorporate the taxa that occur less than the frequency cutoff into
## the "rare" group.
commontaxaT <- taxaT
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


## Remove the list of taxa names that satisfied the frequence cutoff.
rm(freqTaxaT, taxaT)
## ##################################################



## ##################################################
## Write out the data to use in analyses when considering using only
## taxa with at least 5% on at least one subject and day during the
## first 2 weeks.

write.csv(commontaxaT, file="families_hit_cutoff_twice_all_time_steps.csv", row.names=FALSE)
## ##################################################
