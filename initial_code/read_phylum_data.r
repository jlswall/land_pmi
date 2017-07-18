library("readxl")
library("tidyr")
library("tibble")
library("dplyr")
library("stringr")
library("ggplot2")
library("vegan")
library("colorspace")


## library("openxlsx")
## phylumAllT <- read.xlsx(fileNm, sheet="Phylum_all", skipEmptyCols=FALSE)
fileNm <- "../orig_data/Shane_all_skin_samples_taxo_bs_05_05_2017.xlsx"
rawAllT <- read_excel(path=fileNm, sheet="Phylum_all")

## The taxon column contains names of the form "p__name" or "k_name",
## and I will want to make a new column without the "p__" suffix
## later.  So, I rename this column now to avoid confusion later.
colnames(rawAllT)[colnames(rawAllT)=="taxon"] <- "taxonName"

## Some of these column names are duplicated.
## duplicated(colnames(phylumAllDF))
## This first happens at column 116.
## min(which(duplicated(colnames(phylumAllDF))))



## ##################################################
## Extract raw data for each pig in wide format.

## For now, I'm going to completely ignore columns 2 ("rankID") and 4
## ("daughterlevels").  I'm also going to ignore columsn beyond column
## "DJ" (column 114).  These appear to be summary columns.
## The columns I want to focus on are:
## column 3: "taxonName"
## column 5: "total" - to check that the other numbers are correct.
## columns 6-114
##    columns whose names start with "A#" - These are the observed
## counts for the individual pigs at the various times.
##    columns whose names start with "T#" - These are the averages
## of the observed counts over all pigs at the various times.  The
## names of these columns are in form T#days_#accumulatedDegreeDays.

## We're ignoring columns past 114.
mainT <- rawAllT[,1:114]


## Identify column names starting with "A".
namesA <- colnames(mainT)[substring(first=1, last=1, colnames(mainT))=="A"]
wideIndivT <- mainT[,c("taxonName", namesA)]


## Identify column names starting with "T".  The values in these
## columsn are the averages of the individual pigs at each time point.
## The names of these columns contain the number of days since death
## and the accumulated degree days.  The number of days since death
## immediately follows the "T", and the number of accumulated degree
## days follows the "_".
namesT <- colnames(mainT)[substring(first=1, last=1, colnames(mainT))=="T"]
## Separate the days from the accumulated degree days.
timeDF <- separate(data.frame(x=substring(namesT, first=2),
                              stringsAsFactors=F),
                   x, sep="_", into=c("days", "degdays"), convert=T)

## Extract the columns with the phylum names and the average counts
## across pigs for each time point.
wideAvgsT <- mainT[,c("taxonName", namesT)]

rm(namesA, namesT, mainT)
## ##################################################



## ##################################################
## Try to go from wide format to long format.

indivT <- wideIndivT %>%
  gather(indiv_time, counts, -taxonName) %>%
  separate(indiv_time, sep="_T", into=c("subj", "days"), convert=T) %>%
  complete(taxonName, days, subj)
## To see more info about how this works, see:
##   http://www.milanor.net/blog/reshape-data-r-tidyr-vs-reshape2/


## Next thing to check is whether we can get back the values in
## wideAvgsT when we do a summary for indivT.
chkAvgsT <- indivT %>%
  select(taxonName, days, counts) %>%
  group_by(taxonName, days) %>%
  summarize(avgs=mean(counts, na.rm=TRUE)) %>%
  spread(key=days,  value=avgs)
## Match the names to the timeDF frame.
matchNamesV <- na.omit(match(colnames(chkAvgsT), as.character(timeDF$days)))
chkNamesV <- paste(timeDF[matchNamesV,1], timeDF[matchNamesV,2], sep="_")
colnames(chkAvgsT) <- c("taxonName", paste0("T", chkNamesV))
## Order this in the same order as what I read in from the sheet.
reorderChkT <- chkAvgsT[match(wideAvgsT$taxonName, chkAvgsT$taxonName), colnames(wideAvgsT)]

## Compare with the averages I read in from the sheet.
all.equal(wideAvgsT[,1], reorderChkT[,1])  ## Ensure taxonNames in same order.
## Now check the counts, not the taxa names.
apply(wideAvgsT[,-1] - reorderChkT[,-1], 2, summary)
## There is a problem with column T1_27.  As an example, consider the
## row for p__Firmicutes.
subset(indivT, (days==1) & (taxonName=="p__Firmicutes"), "counts")
apply(subset(indivT, (days==1) & (taxonName=="p__Firmicutes"), "counts"), 2, mean)
subset(wideAvgsT, (taxonName=="p__Firmicutes"), "T1_27")
## Look at all the values for "T1_27".
## cbind(reorderChkT[,"T1_27"], wideAvgsT[,"T1_27"], reorderChkT[,"T1_27"]- wideAvgsT[,"T1_27"] )

rm(chkAvgsT, matchNamesV, chkNamesV)
## ##################################################



## ##################################################
## Make other adjustments to the dataset so that it's easier to use.

## Check that all the phylum counts add up to the k__Bacteria counts.
tmp1 <- indivT %>%
  filter(taxonName!="k__Bacteria") %>%
  group_by(days, subj) %>%
  summarize(counts=sum(counts))
tmp2 <- subset(indivT, taxonName=="k__Bacteria", c("days", "subj", "counts"))
if (!all.equal(tmp1, tmp2))
  stop("Phylum counts don't add up to k__Bacteria counts")
rm(tmp1, tmp2)


## Remove the k__Backteria row, since it is just the totals of the
## phyla.  Also, include accum. degree days in the tibble.
indivT <- indivT %>%
  filter(taxonName!="k__Bacteria") %>%
  left_join(timeDF, by="days")


## Make a new, more readable taxa column by stripping the "p__"
## from the beginning of each taxon name.
indivT$taxa <- gsub(indivT$taxonName, pattern="p__", replacement="")
## Remova the taxonName column from the tibble to avoid confusion with
## the next taxa column.
indivT <- indivT %>% select(-taxonName)
## ##################################################



## ##################################################
## For each individual pig and each day, find the percentage of all
## counts represented by each taxa.

ctsByDaySubjT <- indivT %>%
  group_by(days, degdays, subj) %>%
  summarize(totals = sum(counts))

indivT <- indivT %>%
  left_join(ctsByDaySubjT) %>%
  mutate(percByDaySubj=counts/totals) %>%
  select(days, degdays, subj, taxa, counts, percByDaySubj)
## ##################################################



## ##################################################
## Treat each pig as a "site"/"community" and each day as a time step.
## The repsonse is multivariate because there are counts for many
## different phyla.

## We need to get the data into "community matrix format", which is
## basically a wide format.
wideindivT <- indivT %>%
  select(days, degdays, subj, taxa, counts) %>%
  spread(taxa, counts)
## Now, break this into 2 pieces.  One piece is the "community
## matrix", which has a column for each phylum and with the count for
## that phylum for a given day and subject.  The second piece is the
## covariates (subject and day).
communityCovariates <- wideindivT[,1:3]
communityCounts <- wideindivT[,-c(1:3)]
## If the count is NA, make it zero.
communityCounts[is.na(communityCounts)] <- 0
## Get rid of rows which have zero totals.
zeroRows <- apply(communityCounts, 1, sum)==0
communityCounts <- communityCounts[!zeroRows,]
communityCovariates <- communityCovariates[!zeroRows,]
rm(zeroRows)

## Fit the PERMANOVA model
## First without strata:
## trynostrata <- adonis(communityCounts ~ as.factor(days), data=communityCovariates, permutations=5000)
## Then, with strata.  (See
## http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html).
trystrata <- adonis(communityCounts ~ as.factor(days), data=communityCovariates, strata=communityCovariates$subj, permutations=5000)


## Try using nonmetric multidimensional scaling.
trymds <- metaMDS(communityCounts, k=2, trymax=1000)
stressplot(trymds)
plot(trymds)
ordiplot(trymds,type="n")
orditorp(trymds, display="species", col="red", air=0.01,
         labels=make.cepnames(names(communityCounts)))
## Make colors.
myColors <- diverge_hcl(length(unique(communityCovariates$days)))
orditorp(trymds, display="sites", air=0.01, col=myColors[as.numeric(as.factor(communityCovariates$days))], labels=as.character(communityCovariates$days))
## ordihull(trymds, groups=communityCovariates$days, draw="polygon", col=myColors)
rm(myColors)

rm(trystrata, trymds)
## ##################################################




## ##################################################
## Make some exploratory graphics.

## Number of days and accum. degree days are strongly correlated.
with(indivT, cor(degdays, days))

## For each bacteria, plot counts for each pig vs. accum. degree days.
ggplot(indivT, aes(degdays, counts)) +
  geom_point() +
  facet_wrap(~taxa)


## #####################
## Make various stacked bar charts.

## For each day, make stacked bar with layer for each phylum.
## To make this easier to read, any taxa that represent less than 3%
## of the total will be classified as "rare taxa".
rareT <- indivT %>%
    group_by(taxa) %>%
    summarize(taxatotal = sum(counts, na.rm=TRUE)) %>%
    mutate(freq = taxatotal/sum(taxatotal)) %>%
    filter(freq < 0.03) %>%
    select(taxa)
## Rename taxa that occur less than 3% of the time to
## "rare/unclassif. taxa".  The sum the counts for all these
## rare/unclassif. taxa into one row.
renameT <- indivT
renameT[renameT$taxa %in% rareT[[1]], "taxa"] <- "Rare/unclassif. taxa"
renameT <- renameT %>%
  group_by(days, degdays, subj, taxa) %>%
  summarize(counts = sum(counts, na.rm=TRUE))
rm(rareT)

## Now summarize the counts by day and by day-taxa.
barchartT <- renameT %>%
  group_by(days, degdays, taxa) %>%
  summarize(taxadaytotal = sum(counts, na.rm=TRUE)) %>%
  left_join(
      renameT %>%
      group_by(degdays) %>%
      summarize(daytotal = sum(counts, na.rm=TRUE)),
      by = "degdays"
  ) %>%
  mutate(perctaxaday = taxadaytotal / daytotal)

## Make the stacked bar chart using raw counts by degree days.
ggplot(renameT, aes(degdays)) +
  geom_bar(aes(weight=counts, fill=taxa))

## Make stacked bar chart using percentages of taxa for each degree
## day based on the total count for that day across all pigs.
ggplot(barchartT, aes(degdays)) +
  geom_bar(aes(weight=perctaxaday, fill=taxa)) +
  labs(x="Accum. degree days", y="Relative abundance (%)")

## Make stacked bar chart to compare with Figure 1 in Pechal et al
## (2013).
ggplot(subset(barchartT, days <= 5), aes(days)) +
  geom_bar(aes(weight=perctaxaday, fill=taxa)) +
  labs(x="Days", y="Relative abundance (%)")


## Assess variability among pigs on the first few days.
ggplot(subset(renameT, days <= 5), aes(x=subj, y=counts, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days)

## The counts for each day and each subject are really different, so
## we plot these in terms of the proportion per days.
## First find out which taxa aren't considered "rare"; i.e. have more
## than 3% for at least one subject and one day.
commontaxaT <- indivT %>%
  filter(percByDaySubj >= 0.03) %>%
  distinct(taxa)
renameT <- indivT
renameT[!(renameT$taxa %in% commontaxaT[[1]]), "taxa"] <- "Rare taxa"
ggplot(renameT %>%
       group_by(days, degdays, subj, taxa) %>%
       summarize(perc=sum(percByDaySubj)) %>%
       filter(days <= 5),
       aes(x=subj, y=perc, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days)
  
rm(barchartT, renameT)
## #####################


## For each bacteria, plot counts for each pig vs. accum. degree days.
ggplot(subset(indivT, taxa %in% commonTaxa), aes(degdays, counts)) +
  geom_point() +
  facet_wrap(~taxa)
## ##################################################
