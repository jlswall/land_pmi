library("tidyverse")
library("readxl")
library("stringr")


## Read in data from Excel.
fileNm <- "Skin_Validation_Summary_summer_2017.xlsx"
rawAllT <- read_excel(path=fileNm, sheet="Order")
rm(fileNm)


## ##################################################
## Deal with column names.

## The "taxon" column contains some names that are strange because
## they have "k__" or "o__" at the beginning.  I will want to make a
## new column later with only text-based names.  So, I rename this
## column now to avoid confusion later.
colnames(rawAllT)[colnames(rawAllT)=="taxon"] <- "origName"

## Ensure column names are unique in this sheet (unlike the
## "Phylum_all" sheet in the first spreadsheet I got.)
sum(duplicated(colnames(rawAllT)))

## The column names for the individual subject day counts all end with
## "S1" (positions 5-6 in the string).  This makes them harder to deal
## with later in the code, so we remove the "S1".
indicS1 <- substring(colnames(rawAllT), first=5, last=6)=="S1"
colnames(rawAllT)[indicS1] <- substring(colnames(rawAllT)[indicS1], first=1, last=4)
rm(indicS1)
## ##################################################



## ##################################################
## Focus on the counts for the individual pigs on the various
## days.  Want to transfer from wide format to long format.  Check that the
## columns whose names start with "Avg" are actually the averages I
## think they are, and I'll check that totals column and row seems
## correct.

## Column 3: taxa name (some with strange "o__", "c__" prefixes
##   The first row contains the totals.
## Column 5: contains totals of all counts (all days and subjects) for
## the taxa in that row
##
## Columns 6-19: 
## - columns whose names start with "P" - These are the observed
## counts for the individual pigs at the various times.
## - columns whose names start with "Avg_" - These are the averages
## of the observed counts over all pigs at the various times.


## #######################
## Put individual counts and average counts into different tables.

## Identify column names starting with "P". Save these as the counts
## for individual pigs on the various days.
namesP <- colnames(rawAllT)[substring(first=1, last=1, colnames(rawAllT))=="P"]
## Extract the columns with the taxa names and the counts by subject and day.
wideIndivT <- rawAllT[,c("origName", namesP)]

## Identify column names starting with "Avg_T".  The values in these
## columns are the averages of the individual pigs at each time point.
## The names of these columns contain the number of days since death.
## The number of days since death immediately follows the "T".
namesT <- colnames(rawAllT)[substring(first=1, last=5, colnames(rawAllT))=="Avg_T"]
## Extract the columns with the taxa names and the average counts by day.
wideAvgsT <- rawAllT[,c("origName", namesT)]

rm(namesP, namesT)
## #######################


## #######################
## Go from wide format to long format.  Check the averages and totals
## columns.

## Go from wide to long format.  Pull out subject and day from the
## column names of the form "P#T#".
indivT <- wideIndivT %>%
  gather(indiv_time, counts, -origName) %>%
  separate(indiv_time, sep="T", into=c("subj", "days"), convert=T)


## Check whether we can get back the values in wideAvgsT when we do a
## summary for indivT.
chkAvgsT <- indivT %>%
  select(origName, days, counts) %>%
  group_by(origName, days) %>%
  summarize(avgs=mean(counts)) %>%
  spread(key=days,  value=avgs)
## Match the names to the original column names.
matchNamesV <- colnames(chkAvgsT) %in% as.character(unique(indivT$days))
chkNamesV <- paste0("Avg_T", names(chkAvgsT)[matchNamesV])
colnames(chkAvgsT) <- c(colnames(chkAvgsT)[!matchNamesV], chkNamesV)
rm(matchNamesV, chkNamesV)
## Order this in the same order as what I read in from the sheet.
reorderChkT <- chkAvgsT[match(wideAvgsT$origName, chkAvgsT$origName), colnames(wideAvgsT)]


## Compare with the averages I read in from the sheet.
## First, ensure taxonNames in same order.
all.equal(wideAvgsT[,1], reorderChkT[,1])
## Now check the counts.  The averages in the spreadsheet appear to
## have been rounded, so we expect some differences.  However, there
## shouldn't be any differences greater than 0.5.
indicBigDiff <- abs(wideAvgsT[,-1] - reorderChkT[,-1]) > 0.5
sum(indicBigDiff)

rm(indicBigDiff)
## #######################



## #######################
## Check that the total counts for each subject on each day match the
## "Bacteria" row (first row).

## Save these totals in a table for use in calculating percentages.
## We exclude "Bacteria" taxa because that line is supposed to contain
## the totals of all the taxa, including the unclassified taxa.
## Later, I'll re-do these counts to exclude the "unclassified" taxa.
ctBySubjDayT <- indivT %>%
  filter(origName!="k__Bacteria") %>%
  group_by(days, subj) %>%
  summarize(totals=sum(counts))

## Compare the totals calculated above with the first row
## ("k__Bacteria")
all.equal(
    unite(ctBySubjDayT, subj_day, subj, days, sep="T") %>%
      spread(key=subj_day, value=totals),
    wideIndivT %>% filter(origName=="k__Bacteria") %>% select(-origName)
)
## #######################


## #######################
## Check that the total counts for each taxa (across days and pigs)
## match the totals column (column "E").

## There's a PROBLEM HERE, because the totals in column e are often
## bigger than my totals (never less than my totals).
totalDiffT <- indivT %>%
  group_by(origName) %>%
  filter(origName!="k__Bacteria") %>%
  summarize(myTotal=sum(counts)) %>%
  full_join(rawAllT %>% filter(origName!="k__Bacteria") %>% arrange(origName) %>% select(origName, total)) %>%
  mutate(mineMinusTheirs = myTotal-total)

summary(totalDiffT$mineMinusTheirs)
## #######################

rm(wideAvgsT, wideIndivT, reorderChkT, totalDiffT)
## ##################################################



## ##################################################
## Find the percentage of counts which are unclassified bacteria.

## About 0.35% are unclassified.
sum(subset(indivT, origName=="k__Bacteria_unclassified")[,"counts"])/sum(subset(indivT, origName!="k__Bacteria")[,"counts"])
## ##################################################



## ##################################################
## Make other adjustments to the dataset so that it's easier to use.

## Read in the table containing degdays by subj and time.
timeT <- read_csv("degdays_by_subj_day.csv")


## Remove the k__Bacteria row (first row), since it is just the totals
## of the taxa.  Remove the counts associated with unclassifed taxa.
## Also, include accum. degree days in the tibble.
indivT <- indivT %>%
  filter(!(origName %in% c("k__Bacteria", "k__Bacteria_unclassified"))) %>%
  left_join(timeT)


## ## Make a new, more readable taxa column.
## ## Column names with open brackets (e.g. "[Tissierellaceae]") causes
## ## problems for functions expecting traditional data frame column
## ## names.
## indivT$taxa <- gsub(indivT$origName, pattern="\\[", replacement="")
## indivT$taxa <- gsub(indivT$taxa, pattern="]", replacement="")
## ## Column names with dashes can likewise be a problem, so I replace
## ## dashes with underscores.
## indivT$taxa <- gsub(indivT$taxa, pattern="-", replacement="_")

## I decided to leave the names alone for the taxa, even though they
## aren't very read-able.
indivT$taxa <- indivT$origName
## Remove the taxonName column from the tibble to avoid confusion with
## the next taxa column.
indivT <- indivT %>% select(-origName)
## ##################################################



## ##################################################
## For use in graphs and in calculating percentages later, we need
## total counts (over all taxa, unclassified taxa excluded) by:
##   Each pig and each day 
##   Each day (all pigs combined)

## Total taxa counts by day and subject (each pig separately).
ctBySubjDayT <- indivT %>%
  group_by(days, degdays, subj) %>%
  summarize(totals=sum(counts))

## Total taxa counts by day (all pigs combined).
ctByDayT <- indivT %>%
  group_by(days, degdays) %>%
  summarize(totals = sum(counts))
## ##################################################



## ##################################################
## Some taxa don't occur frequently.  It's hard to make a hard cutoff
## for what constitutes "frequently".  There are 140 taxa in the
## dataset, and a lot of them appear in less than 0.1% of samples.

## I'm going to set the cutoff at 1% (0.01).  This means that in order
## to be included in the dataset, a specific taxa must make up at
## least 1% of the total counts on at least 1 day for at least 1
## cadaver.
freqCutoff <- 0.01

## Get list of maximum taxa percentages sorted in descending order:
data.frame(indivT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay = counts/totals) %>%
  group_by(taxa) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  arrange(desc(maxFracBySubjDay))
)


## Save the taxa names (in a tibble) which satisfy the frequency
## cutoff.
freqTaxaT <- indivT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay = counts/totals) %>%
  group_by(taxa) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  filter(maxFracBySubjDay >= freqCutoff) %>%
  arrange(desc(maxFracBySubjDay)) %>%
  select(taxa)


## Rename taxa that occur less than the frequency cutoff allows as
## "rare".  Then, sum all these "rare" taxa into one row.
commontaxaT <- indivT
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



## ############ WORKING HERE! ############



## ##################################################
## Save the tibble to a file for use in separate code
## for graphing and analysis.

## I have to use the base R write.csv() routine, because write_csv
## will write out scientific notation, which read_csv() doesn't read
## in properly.
## write_csv(commontaxaT, path="orders_massaged.csv")
write.csv(commontaxaT, file="orders_massaged.csv", row.names=FALSE)
## ##################################################
