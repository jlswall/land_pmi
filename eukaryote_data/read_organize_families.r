library("tidyverse")
library("readxl")
library("stringr")


## Read in data from Excel.
fileNm <- "luisa_updated_2018-04-29.xlsx"
rawAllT <- read_excel(path=fileNm, sheet="Rank 5 - Family", skip=1)
rm(fileNm)

## Column names are unique in this sheet.
sum(duplicated(colnames(rawAllT)))



## ##################################################
## Concentrate on the counts for the individual pigs on the various
## days.  Transfer from wide format to long format.  Check that the
## columns whose names start with "T" are actually the sums I
## think they are, and I'll check that totals column and row seems
## correct.

## #######################
## Put individual counts and daily sums into different tables.

## Identify column names starting with "A". Save these as the counts
## for individual pigs on the various days.
namesA <- colnames(rawAllT)[substring(first=1, last=1, colnames(rawAllT))=="A"]
wideIndivT <- rawAllT[,c("taxon", namesA)]

## Identify column names starting with "T".  The values in these
## columns are the sums of the individual pigs at each time point.
## The names of these columns contain the number of days since death
## and the accumulated degree days.  The number of days since death
## immediately follows the "T", and the number of accumulated degree
## days follows the "_".
namesT <- colnames(rawAllT)[substring(first=1, last=1, colnames(rawAllT))=="T"]
## Separate the days from the accumulated degree days.
timeDF <- separate(data.frame(x=substring(namesT, first=2),
                              stringsAsFactors=F),
                   x, sep=" - ", into=c("days", "degdays"), convert=T)
## Note: number of days and accum. degree days are strongly correlated.
with(timeDF, cor(degdays, days))


## Extract the columns with the taxa names and the sums across pigs
## for each time point and taxon.
wideSumsT <- rawAllT[,c("taxon", namesT)]


## The last few rows are special cases, so I exclude them from the
## tables of counts for the taxa.  They are: "Eukaryota",
## '"Unclassified"', "% Unclassified", "Total Classified"
wideIndivT <- wideIndivT %>% filter(!(taxon %in% c('Eukaryota', '"Unclassified"', "% Unclassified", "Total Classified")))
wideSumsT <- wideSumsT %>% filter(!(taxon %in% c('Eukaryota', '"Unclassified"', "% Unclassified", "Total Classified")))


rm(namesA, namesT)
## #######################


## #######################
## Go from wide format to long format.  Check the columns and rows
## containing totals.

## Go from wide to long format, separating column names into subject,
## day, and swab number.
indivT <- wideIndivT %>%
  gather(indiv_time, counts, -taxon) %>%
  separate(indiv_time, sep="T", into=c("subj", "days"), convert=T) %>%
  separate(days, sep="S", into=c("days", "swab"), convert=T)


## ######
## Make sure the totals for each day and taxa (across subjects) match
## those I get from doing the sums.
mysumsT <- indivT %>% group_by(taxon, days) %>% summarize(counts=sum(counts))
## Adjust the times so that they match up with the first few
## characters of the original column names for comparison.
mysumsT[,"compareDays"] <- paste0("T", as.vector(as.matrix(mysumsT[,"days"])), " - ")
mysumsT <- mysumsT %>% select(-days) %>% spread(compareDays, counts)
## Put column names of mysumsT into same order as those for
## wideSums. It's sufficient to just check that the "T##" part
## matches.
match.order <- match(substring(colnames(wideSumsT), first=1, last=4), substring(colnames(mysumsT), first=1, last=4))
mysumsT <- mysumsT[,match.order]
## Put rows of mysumsT into same order as those for wideSums.
match.order <- match(wideSumsT$taxon, mysumsT$taxon)
mysumsT <- mysumsT[match.order,]
unique(as.vector(as.matrix(wideSumsT[,-1]) - as.matrix(mysumsT[,-1])))

rm(match.order, wideSumsT, mysumsT)
## ######
## #######################
## ##################################################




## ##################################################
## Find the total percentage of counts which are unclassified.

## About 32.76% are unclassified.
pull(indivT %>%
     filter(grepl("_unclassified|_uncultured|Incertae_Sedis", taxon)) %>%
     summarize(total_uncl=sum(counts)), "total_uncl") / sum(indivT[,"counts"])
## ##################################################




## ##################################################
## Make other adjustments to the dataset so that it's easier to use.

## Remove the counts associated with unclassifed taxa,
## "Incertae_Sedis", and "uncultured" (*_unclassified, "_uncultured",
## "Incertae_Sedis").  Also, include accum. degree days in the tibble.
indivT <- indivT %>%
  filter(!grepl("_unclassified|_uncultured|Incertae_Sedis", taxon)) %>%
  left_join(timeDF, by="days")

## Remove the counts associated with order "Mammalia" (which is
## probably the pig's DNA), and with order "Aves" (birds, not sure why
## that's in there).
indivT <- indivT %>%
  filter(!grepl("Mammalia|Aves", taxon))
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
## for what constitutes "frequently".  There are 164 (classified) taxa
## in the dataset, and a lot of them appear in less than 0.1% of
## samples.

## I'm going to set the cutoff at 1% (0.01).  This means that in order
## to be included in the dataset, a specific taxa must make up at
## least 1% of the total counts on at least 1 day for at least 1
## cadaver.
freqCutoff <- 0.01

## Get list of maximum taxa percentages sorted in descending order:
data.frame(indivT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay = counts/totals) %>%
  group_by(taxon) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  arrange(desc(maxFracBySubjDay))
)


## Save the taxa names (in a tibble) which satisfy the frequency
## cutoff.
freqTaxaT <- indivT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay = counts/totals) %>%
  group_by(taxon) %>%
  summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
  filter(maxFracBySubjDay >= freqCutoff) %>%
  arrange(desc(maxFracBySubjDay)) %>%
  select(taxon)


## Rename taxa that occur less than the frequency cutoff allows as
## "rare".  Then, sum all these "rare" taxa into one row.
commontaxaT <- indivT
commontaxaT[!(commontaxaT$taxon %in% freqTaxaT$taxon), "taxon"] <- "Rare"
commontaxaT <- commontaxaT %>%
  group_by(days, degdays, subj, taxon) %>%
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
write.csv(commontaxaT, file="families_massaged.csv", row.names=FALSE)
## ##################################################
