library("tidyverse")

## ##################################################
## Read the eukaryote and bacteria data which we previously organized
## separately.
eukT <- read_csv("../eukaryote_data/orders_massaged.csv")
bacT <- read_csv("../revise_algorithm/orders_massaged.csv")

## Combine these two tibbles, leaving out the column with the
## fractional calculations and renaming the "taxon" column in the
## eukaryote data to "taxa".
comboT <- rbind(eukT %>% rename(taxa=taxon) %>% select(-fracBySubjDay),
                bacT %>% select(-fracBySubjDay))


## For the bacteria data, we did not have data for
## - subject A1 on days 7 and 9
## - subject A4 on day 7
## We did have data for these subjects and days in the eukaryote data,
## but because of the missing bacteria data, we remove these rows from
## the tibble.
comboT <- comboT %>% filter(! ((subj=="A1") & (days==7 | days==9)) )
comboT <- comboT %>% filter( !((subj=="A4") & (days==7)) )

## For the bacteria data, we had problems with the data from
## - subject A3 on day 40
## We did have data for these subjects and days in the eukaryote data,
## but because of the missing bacteria data, we remove these rows from
## the tibble.
comboT <- comboT %>% filter( !((subj=="A3") & (days==40)) )
## ##################################################



## ##################################################
## Total taxa counts by day and subject (each pig separately).
ctBySubjDayT <- comboT %>%
  group_by(days, degdays, subj) %>%
  summarize(totals=sum(counts))
## ##################################################



## ##################################################
## Some taxa don't occur frequently.  It's hard to make a hard cutoff
## for what constitutes "frequently".  However, if we leave all 230
## taxa in the dataset, things become overwhelming with taxa that make
## up, for example, less than even 0.1% of the counts for a specific
## cadaver on a specific day.

## I'm going to set the cutoff at 1% (0.01).  This means that in order
## to be included in the dataset, a specific taxa must make up at
## least 1% of the total counts on at least 1 day for at least 1
## cadaver.
freqCutoff <- 0.01

## ## Get list of maximum taxa percentages sorted in descending order:
## data.frame(comboT %>%
##   left_join(ctBySubjDayT) %>%
##   mutate(fracBySubjDay = counts/totals) %>%
##   group_by(taxa) %>%
##   summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
##   arrange(desc(maxFracBySubjDay))
## )

## Save the taxa names (in a tibble) which satisfy the frequency
## cutoff.
freqTaxaT <- comboT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay = counts/totals) %>%
  group_by(taxa) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  filter(maxFracBySubjDay >= freqCutoff) %>%
  arrange(desc(maxFracBySubjDay)) %>%
  select(taxa)


## Rename taxa that occur less than the frequency cutoff allows as
## "rare".  Then, sum all these "rare" taxa into one row.
commontaxaT <- comboT
commontaxaT[!(commontaxaT$taxa %in% freqTaxaT$taxa), "taxa"] <- "Rare"
commontaxaT <- commontaxaT %>%
  group_by(days, degdays, subj, taxa) %>%
  summarize(counts = sum(counts))

## Remove the list of taxa names that satisfied the frequence cutoff.
rm(freqTaxaT)
## ##################################################



## ##################################################
## Add percentages by subj/day to the commontaxaT table.

## Use the table of total counts by subj/day to find the fraction
## represented by each taxa for each subj/day.
commontaxaT <- commontaxaT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay=counts/totals) %>%
  select(-totals)

## Check that the fractions add up to 1, appropriately.
unique(
    unlist(commontaxaT %>%
           group_by(days, subj) %>%
           summarize(sumFracBySubjDay = sum(fracBySubjDay)) %>%
           ungroup() %>%
           select(sumFracBySubjDay))
    )
## ##################################################



## ##################################################
## Save the tibble to a file for use in separate code
## for graphing and analysis.

## I decided to use the write.csv() routine for the order data, even
## though I didn't have the problems here that I did with the orders
## and the phyla.  That way, it remains consistent, and hopefully I
## prevent future issues with write_csv and scientific notation.
write.csv(commontaxaT, file="orders_massaged.csv", row.names=FALSE)
## ##################################################
