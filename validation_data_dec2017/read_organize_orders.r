library("tidyverse")
library("readxl")
library("stringr")


## Read in data from Excel.
fileNm <- "Skin_Validation_Summary_summer_2017.xlsx"
rawAllT <- read_excel(path=fileNm, sheet="Order")
rm(fileNm)

## The "taxon" column contains some names that are strange because
## they have "k__" or "o__" at the beginning.  I will want to make a
## new column later with only text-based names.  So, I rename this
## column now to avoid confusion later.
colnames(rawAllT)[colnames(rawAllT)=="taxon"] <- "origName"

## Column names are unique in this sheet (unlike the "Phylum_all"
## sheet in the original data I got.)
sum(duplicated(colnames(rawAllT)))



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
## column names of the form "P#T#S1".
indivT <- wideIndivT %>%
  gather(indiv_time, counts, -origName) %>%
  separate(indiv_time, sep="T", into=c("subj", "days"), convert=T) %>%
  separate(days, sep="S", into=c("days", "notneeded"), convert=T) %>%
  select(-notneeded)


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


## ########## THERE ARE SOME MIS-MATCHES BELOW ##########
## ########## I think the averages in the sheet are rounded.

## Compare with the averages I read in from the sheet.
## First, ensure taxonNames in same order.
all.equal(wideAvgsT[,1], reorderChkT[,1])
## Now check the counts, not the taxa names.
apply(wideAvgsT[,-1] - reorderChkT[,-1], 2, summary)
## There is a problem with column T1_27.  As an example, consider the
## counts for Clostridiales on that day.
subset(indivT, (days==1) & (origName=="Clostridiales"), "counts")
## The average is given by:
apply(subset(indivT, (days==1) & (origName=="Clostridiales"), "counts"), 2, mean)
## But, this is not the average we read from the sheet.
subset(wideAvgsT, (origName=="Clostridiales"), "T1_27")
## It looks like they took the sum, not the mean.
apply(subset(indivT, (days==1) & (origName=="Clostridiales"), "counts"), 2, sum)
## Look at all the values for "T1_27".
## cbind(reorderChkT[,"T1_27"], wideAvgsT[,"T1_27"], reorderChkT[,"T1_27"]- wideAvgsT[,"T1_27"] )
rm(chkAvgsT, matchNamesV, chkNamesV)
## #######################


## #######################
## Check that the total counts for each taxa match the "total" column
## (column #111).  There are a lot of these to check, so we take the
## absolute differences between the our total counts and the "total"
## column and make sure that the biggest difference is 0.
max(indivT %>%
  group_by(origName) %>%
  summarize(totalCt = sum(counts)) %>%
  left_join(rawAllT %>% select(origName, total)) %>%
  mutate(absDiffOrigMyCalc = abs(total - totalCt)) %>%
  select(absDiffOrigMyCalc)
  )
## #######################


## #######################
## Check that the total counts for each subject on each day match the
## "Bacteria" row (row #233 in the tibble, #234 in the worksheet).

## Save these totals in a table for use in calculating percentages.
## We exclude "Bacteria" taxa because that line is supposed to contain
## the totals of all the taxa, including the unclassified taxa.
## Later, I'll re-do these counts to exclude the "unclassified" taxa.
ctBySubjDayT <- indivT %>%
  filter(origName!="Bacteria") %>%
  group_by(days, subj) %>%
  summarize(totals=sum(counts))

## Compare the totals calculated above with the last row ("Bacteria")
## of the individual counts.
all.equal(
    unite(ctBySubjDayT, subj_day, subj, days, sep="_T") %>%
      spread(key=subj_day, value=totals),
    wideIndivT %>% filter(origName=="Bacteria") %>% select(-origName)
)
## #######################

rm(mainT, wideAvgsT, wideIndivT, reorderChkT)
## ##################################################




## ##################################################
## Columns 113-223 ("DI"-"HO") appear to be percentages for each taxa, by
## day and pig.  Check these.
## Columns 225-242: I haven't checked these yet.

## #######################
## Organize the information.

## Put these columns in their own table.
widePercT <- rawAllT[,113:223]

## The first column contains the order names.
colnames(widePercT)[1] <- "origName"

## Identify column names starting with "A" (individuals A1-A6).
namesA <- colnames(widePercT)[substring(first=1, last=1, colnames(widePercT))=="A"]
wideIndivPercT <- rawAllT[,c("origName", namesA)]


## Identify column names starting with "T" (averages across
## individuals).
namesT <- colnames(widePercT)[substring(first=1, last=1, colnames(widePercT))=="T"]
wideAvgsPercT <- rawAllT[,c("origName", namesA)]

rm(namesA, namesT)
## #######################


## #######################
## Try to take the individual percentages from wide format to long
## format.
## indivPercT <- wideIndivPercT %>%
##   gather(indiv_time, perc, -origName) %>%
##   separate(indiv_time, sep="_T", into=c("subj", "days_with_extra"), convert=T) %>%
##   separate(days_with_extra, sep="__", into=c("days", "extra_stuff"), convert=T) %>%
##   select(-extra_stuff)
## #######################


## #######################
## Check that these percentages straight from the worksheet are the
## same as what we calculate based on the individual counts.

## First, I calculate these percentages based on the individuals
## counts and the sums (based on those counts) that I calculated
## earlier.  I add this column to the main table.
indivT <- indivT %>%
  left_join(ctBySubjDayT) %>%
  mutate(percByDaySubj = 100*counts/totals) %>%
  select(-totals)


## Try to put these percentages in wide format for comparison with the
## raw numbers from the worksheet.
chkPercT <- indivT %>%
  select(-counts) %>%
  mutate(extrachar="1") %>%
  unite(subj_days, subj, days, sep="_T") %>%
  unite(subj_days, subj_days, extrachar, sep="__") %>%
  spread(key=subj_days, value=percByDaySubj)

## Now check to see if this matches the numbers we got from the spreadsheet.
if (!all.equal(chkPercT, wideIndivPercT))
  stop("Something different between the two sets of percentages.")
if (nrow(setdiff(chkPercT, wideIndivPercT)) != 0)
  stop("Extra observations were created when working with indivT")
if (nrow(setdiff(wideIndivPercT, chkPercT)) != 0)
  stop("More observations were are in the original worksheet than created when working with indivT")


rm(chkPercT, ctBySubjDayT, wideIndivPercT, widePercT, wideAvgsPercT)
## #######################
## ##################################################




## ##################################################
## Find the percentage of counts which are unclassified.

## About 1.1% are unclassified.
sum(subset(indivT, origName=="Unclassified")[,"counts"])/sum(subset(indivT, origName!="Bacteria")[,"counts"])
## ##################################################




## ##################################################
## Make other adjustments to the dataset so that it's easier to use.

## Remove the Bacteria row (last row), since it is just the totals of
## the taxa.  Remove the counts associated with unclassifed taxa.
## Also, include accum. degree days in the tibble.
indivT <- indivT %>%
  filter(!(origName %in% c("Bacteria", "Unclassified"))) %>%
  left_join(timeDF, by="days")


## Make a new, more readable taxa column.
## Column names with open brackets (e.g. "[Tissierellaceae]") causes
## problems for functions expecting traditional data frame column
## names.
indivT$taxa <- gsub(indivT$origName, pattern="\\[", replacement="")
indivT$taxa <- gsub(indivT$taxa, pattern="]", replacement="")
## Column names with dashes can likewise be a problem, so I replace
## dashes with underscores.
indivT$taxa <- gsub(indivT$taxa, pattern="-", replacement="_")
## Remova the taxonName column from the tibble to avoid confusion with
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
## for what constitutes "frequently".  There are 143 taxa in the
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



## ##################################################
## Save the tibble to a file for use in separate code
## for graphing and analysis.

## I have to use the base R write.csv() routine, because write_csv
## will write out scientific notation, which read_csv() doesn't read
## in properly.
## write_csv(commontaxaT, path="orders_massaged.csv")
write.csv(commontaxaT, file="orders_massaged.csv", row.names=FALSE)
## ##################################################
