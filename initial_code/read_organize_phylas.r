library("tidyverse")
library("readxl")
library("stringr")


## Read in data from the phylum worksheet.
fileNm <- "../orig_data/Shane_all_skin_samples_taxo_bs_05_05_2017.xlsx"
rawAllT <- read_excel(path=fileNm, sheet="Phylum_all")
rm(fileNm)

## The taxon column contains names of the form "p__name" or "k_name",
## and I will want to make a new column without the "p__" suffix
## later.  So, I rename this column now to avoid confusion later.
colnames(rawAllT)[colnames(rawAllT)=="taxon"] <- "origName"


## ##################################################
## Extract raw data for each pig in wide format.

## For now, I'm going to completely ignore columns 2 ("rankID") and 4
## ("daughterlevels").  I'm also going to ignore columsn beyond column
## "DJ" (column 114).  These appear to be summary columns.
## The columns I want to focus on are:
## column 3: "origName"
## column 5: "total" - to check that the other numbers are correct.
## columns 6-114
##    columns whose names start with "A#" - These are the observed
## counts for the individual pigs at the various times.
##    columns whose names start with "T#" - These are the averages
## of the observed counts over all pigs at the various times.  The
## names of these columns are in form T#days_#accumulatedDegreeDays.

## The raw data appears to be in the first 114 columns.
mainT <- rawAllT[,1:114]


## Identify column names starting with "A".
namesA <- colnames(mainT)[substring(first=1, last=1, colnames(mainT))=="A"]
wideIndivT <- mainT[,c("origName", namesA)]


## Identify column names starting with "T".  The values in these
## columns are the averages of the individual pigs at each time point.
## The names of these columns contain the number of days since death
## and the accumulated degree days.  The number of days since death
## immediately follows the "T", and the number of accumulated degree
## days follows the "_".
namesT <- colnames(mainT)[substring(first=1, last=1, colnames(mainT))=="T"]
## Separate the days from the accumulated degree days.
timeDF <- separate(data.frame(x=substring(namesT, first=2),
                              stringsAsFactors=F),
                   x, sep="_", into=c("days", "degdays"), convert=T)
## Note: number of days and accum. degree days are strongly correlated.
with(timeDF, cor(degdays, days))



## Extract the columns with the phylum names and the average counts
## across pigs for each time point.
wideAvgsT <- mainT[,c("origName", namesT)]

rm(namesA, namesT, mainT)
## ##################################################




## ##################################################
## Try to go from wide format to long format.

indivT <- wideIndivT %>%
  gather(indiv_time, counts, -origName) %>%
  separate(indiv_time, sep="_T", into=c("subj", "days"), convert=T)
  ## To add rows of missing values for combinations of days and
  ## subjects on which no samples are available (some subjects
  ## were not observed on certain days), add %>% and then this:
  ## complete(origName, days, subj)


## Next thing to check is whether we can get back the values in
## wideAvgsT when we do a summary for indivT.
chkAvgsT <- indivT %>%
  select(origName, days, counts) %>%
  group_by(origName, days) %>%
  summarize(avgs=mean(counts, na.rm=TRUE)) %>%
  spread(key=days,  value=avgs)
## Match the names to the timeDF frame.
matchNamesV <- na.omit(match(colnames(chkAvgsT), as.character(timeDF$days)))
chkNamesV <- paste(timeDF[matchNamesV,1], timeDF[matchNamesV,2], sep="_")
colnames(chkAvgsT) <- c("origName", paste0("T", chkNamesV))
## Order this in the same order as what I read in from the sheet.
reorderChkT <- chkAvgsT[match(wideAvgsT$origName, chkAvgsT$origName), colnames(wideAvgsT)]

## Compare with the averages I read in from the sheet.
all.equal(wideAvgsT[,1], reorderChkT[,1])  ## Ensure origNames in same order.
## Now check the counts, not the taxa names.
apply(wideAvgsT[,-1] - reorderChkT[,-1], 2, summary)
## There is a problem with column T1_27.  As an example, consider the
## row for p__Firmicutes.
subset(indivT, (days==1) & (origName=="p__Firmicutes"), "counts")
apply(subset(indivT, (days==1) & (origName=="p__Firmicutes"), "counts"), 2, mean)
subset(wideAvgsT, (origName=="p__Firmicutes"), "T1_27")
## It looks like they took the sum, not the mean.
apply(subset(indivT, (days==1) & (origName=="p__Firmicutes"), "counts"), 2, sum)
subset(wideAvgsT, (origName=="p__Firmicutes"), "T1_27")
## Look at all the values for "T1_27".
## cbind(reorderChkT[,"T1_27"], wideAvgsT[,"T1_27"], reorderChkT[,"T1_27"]- wideAvgsT[,"T1_27"] )

rm(chkAvgsT, matchNamesV, chkNamesV)



## Check that all the phylum counts add up to the k__Bacteria counts.
tmp1 <- indivT %>%
  filter(origName!="k__Bacteria") %>%
  group_by(days, subj) %>%
  summarize(counts=sum(counts))
tmp2 <- subset(indivT, origName=="k__Bacteria", c("days", "subj", "counts"))
if (!all.equal(tmp1, tmp2))
  stop("Phylum counts don't add up to k__Bacteria counts")
rm(tmp1, tmp2)
## ##################################################



## ##################################################
## Make other adjustments to the dataset so that it's easier to use.

## Remove the k__Backteria row, since it is just the totals of the
## phyla.  Remove the counts associated with unclassified taxa.  Also,
## include accum. degree days in the tibble.
indivT <- indivT %>%
  filter(!(origName %in% c("k__Bacteria", "unclassified"))) %>%
  left_join(timeDF, by="days")

## Make a new, more readable taxa column by stripping the "p__"
## from the beginning of each taxon name.
indivT$taxa <- gsub(indivT$origName, pattern="p__", replacement="")
## The column name "[Thermi]" causes problems for functions expecting
## traditional data frame column names.
indivT$taxa <- gsub(indivT$taxa, pattern="\\[", replacement="")
indivT$taxa <- gsub(indivT$taxa, pattern="]", replacement="")
## Remova the origName column from the tibble to avoid confusion with
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
  summarize(totals = sum(counts))

## Total taxa counts by day (all pigs combined).
ctByDayT <- indivT %>%
  group_by(days, degdays) %>%
  summarize(totals = sum(counts))
## ##################################################




## ##################################################
## Some taxa don't occur frequently.  We identify these, classify them
## as rare, and then reformat the data accordingly.
## Pechal et al (2013) call "rare" any taxa with <3% of total
## abundance.

## Find taxa which make up >= 3% of the total counts, over all days
## and subjects.  We'll cause these "common" taxa.  Remember that the
## "unclassified" taxa are excluded before we get to this point.
commonByTotalV <- unlist(indivT %>%
    group_by(taxa) %>%
    summarize(taxatotal = sum(counts)) %>%
    mutate(frac = taxatotal/sum(taxatotal)) %>%
    filter(frac >= 0.03) %>%
    select(taxa)
    )
## See a list of all taxa percentages sorted in descending order:
## data.frame( indivT %>%
##     group_by(taxa) %>%
##     summarize(taxatotal = sum(counts)) %>%
##     mutate(frac = taxatotal/sum(taxatotal)) %>%
##     arrange(desc(frac))
##     )


## Taxa counts can vary widely by day and individual.  So, another way
## to figure out which taxa are "common" would be to include any taxa
## making up at least 3% of the total count on at least one particular
## day for any particular subject.  Remember that "unclassified" taxa
## are excluded before we get to this point.
commonBySubjDayV <- unlist(indivT %>%
                           left_join(ctBySubjDayT) %>%
                           mutate(fracBySubjDay = counts/totals) %>%
                           group_by(taxa) %>%
                           summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
                           filter(maxFracBySubjDay >= 0.03) %>%
                           select(taxa)
                           )
## See a list of maximum taxa percentages sorted in descending order:
## data.frame(indivT %>%
##            left_join(ctBySubjDayT) %>%
##            mutate(fracBySubjDay = counts/totals) %>%
##            group_by(taxa) %>%
##            summarize(maxFracBySubjDay = max(fracBySubjDay)) %>%
##            arrange(desc(maxFracBySubjDay))
##            )


## Yet another way to calculate the percentages represented by each
## taxa is to calculate the percentages by day (over all subjects).
## Again, "unclassified" taxa are already excluded by this point.
commonByDayV <- unlist(indivT %>%
                       group_by(days, taxa) %>%
                       summarize(ctByDayTaxa = sum(counts)) %>%
                       left_join(ctByDayT) %>%
                       mutate(fracByDay = ctByDayTaxa/totals) %>%
                       group_by(taxa) %>%
                       summarize(maxFracByDay = max(fracByDay)) %>%
                       filter(maxFracByDay >= 0.03) %>%
                       select(taxa)
                       )
## See a list of maximum taxa by day percentages sorted in descending
## order:
## data.frame(indivT %>%
##   group_by(days, taxa) %>%
##   summarize(ctByDayTaxa = sum(counts)) %>%
##   left_join(ctByDayT) %>%
##   mutate(fracByDay = ctByDayTaxa/totals) %>%
##   group_by(taxa) %>%
##   summarize(maxFracByDay = max(fracByDay)) %>%
##   arrange(desc(maxFracByDay))
##   )


## #######################
## Now, I'm going to look again at the percentages represented by each
## taxa by day/subject.  Then, for each day, I'm going to figure out
## the average percentage represented by each taxa (averaged over all
## subjects).

## First, calculate percentages.
percSubjDayT <- indivT %>%
  left_join(ctBySubjDayT) %>%
  mutate(fracBySubjDay = counts/totals) %>%
  select(-totals)

## For each day, find the mean fraction represented by each taxa
## (averaged over subjects).
avgSubjDayT <- percSubjDayT %>%
  group_by(days, taxa) %>%
  summarize(avgFracByDay = mean(fracBySubjDay),
         sdAvgFracByDay = sd(fracBySubjDay)
         )

## To be classified as a "common" (or not rare taxa), a taxa has to
## average at least 0.03 (3%) for at least one day.
commonByAvgByDayV <-  unlist(
    avgSubjDayT %>%
    filter(avgFracByDay >= 0.03) %>%
    ungroup() %>%
    select(taxa) %>%
    distinct(taxa)
)

rm(percSubjDayT)
## #######################
## ##################################################



## ##################################################
## Figure out which taxa should be treated as "rare".  Re-calculate
## the various counts with all the rare individual taxa grouped into a
## taxa called "Rare".

## NOTE: For phylum, lists of common family-level taxa are the same,
## regardless of how we calculate them.
all.equal(commonByAvgByDayV, commonByDayV, commonBySubjDayV, commonByTotalV)
commonTaxaNamesV <- commonByAvgByDayV

rm(commonBySubjDayV, commonByDayV, commonByTotalV, commonByAvgByDayV)
rm(avgSubjDayT)


## Rename taxa that occur less than 3% of the time to "rare".  Then,
## sum the counts for all these rare or uncl. taxa into one row.
commontaxaT <- indivT
commontaxaT[!(commontaxaT$taxa %in% commonTaxaNamesV), "taxa"] <- "Rare"
commontaxaT <- commontaxaT %>%
  group_by(days, degdays, subj, taxa) %>%
  summarize(counts = sum(counts))
rm(commonTaxaNamesV)
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

write_csv(commontaxaT, path="phyla_massaged.csv")
## ##################################################
