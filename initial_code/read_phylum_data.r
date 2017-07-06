library("readxl")
library("tidyr")
library("tibble")
library("dplyr")
library("stringr")
library("ggplot2")


## library("openxlsx")
## phylumAllT <- read.xlsx(fileNm, sheet="Phylum_all", skipEmptyCols=FALSE)
fileNm <- "../orig_data/Shane_all_skin_samples_taxo_bs_05_05_2017.xlsx"
rawAllT <- read_excel(path=fileNm, sheet="Phylum_all")

## Some of these column names are duplicated.
## duplicated(colnames(phylumAllDF))
## This first happens at column 116.
## min(which(duplicated(colnames(phylumAllDF))))


## The first column contains "2"s for all the rows, except the very
## last row (row #36).
unique(rawAllT[1:35,1])
rawAllT[36, 1]
## The last row appears to be at the kingdom level, while the
## preceding rows were at the phylum level.


## ##################################################
## Extract raw data for each pig in wide format.

## For now, I'm going to completely ignore columns 2 ("rankID") and 4
## ("daughterlevels").  I'm also going to ignore columsn beyond column
## "DJ" (column 114).  These appear to be summary columns.
## The columns I want to focus on are:
## column 3: "taxon"
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
wideIndivT <- mainT[,c("taxon", namesA)]


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
wideAvgsT <- mainT[,c("taxon", namesT)]

rm(namesA, namesT, mainT)
## ##################################################



## ##################################################
## Try to go from wide format to long format.

indivT <- wideIndivT %>%
  gather(indiv_time, counts, -taxon) %>%
  separate(indiv_time, sep="_T", into=c("subj", "days"), convert=T) %>%
  complete(taxon, days, subj)
## To see more info about how this works, see:
##   http://www.milanor.net/blog/reshape-data-r-tidyr-vs-reshape2/


## Next thing to check is whether we can get back the values in
## wideAvgsT when we do a summary for indivT.
chkAvgsT <- indivT %>%
  select(taxon, days, counts) %>%
  group_by(taxon, days) %>%
  summarize(avgs=mean(counts, na.rm=TRUE)) %>%
  spread(key=days,  value=avgs)
## Match the names to the timeDF frame.
matchNamesV <- na.omit(match(colnames(chkAvgsT), as.character(timeDF$days)))
chkNamesV <- paste(timeDF[matchNamesV,1], timeDF[matchNamesV,2], sep="_")
colnames(chkAvgsT) <- c("taxon", paste0("T", chkNamesV))
## Order this in the same order as what I read in from the sheet.
reorderChkT <- chkAvgsT[match(wideAvgsT$taxon, chkAvgsT$taxon), colnames(wideAvgsT)]

## Compare with the averages I read in from the sheet.
all.equal(wideAvgsT[,1], reorderChkT[,1])  ## Ensure taxons in same order.
## Now check the counts, not the taxa names.
apply(wideAvgsT[,-1] - reorderChkT[,-1], 2, summary)
## There is a problem with column T1_27.  As an example, consider the
## row for p__Firmicutes.
subset(indivT, (days==1) & (taxon=="p__Firmicutes"), "counts")
apply(subset(indivT, (days==1) & (taxon=="p__Firmicutes"), "counts"), 2, mean)
subset(wideAvgsT, (taxon=="p__Firmicutes"), "T1_27")
## Look at all the values for "T1_27".
cbind(reorderChkT[,"T1_27"], wideAvgsT[,"T1_27"], reorderChkT[,"T1_27"]- wideAvgsT[,"T1_27"] )

rm(chkAvgsT, matchNamesV, chkNamesV)
## ##################################################



## ##################################################
## Check that all the phylum counts add up to the k__Bacteria counts.

tmp1 <- indivT %>%
  filter(taxon!="k__Bacteria") %>%
  group_by(days, subj) %>%
  summarize(counts=sum(counts))

tmp2 <- subset(indivT, taxon=="k__Bacteria", c("days", "subj", "counts"))

all.equal(tmp1, tmp2)

rm(tmp1, tmp2)
## ##################################################



## ##################################################
## Make some exploratory graphics.

## Include accum. degree days in the tibble.
indivT <- indivT %>%
  filter(taxon!="k__Bacteria") %>%
  left_join(timeDF, by="days")

## Make a new species column by stripping the "p__" from the beginning
## of each taxon name.
indivT$species <- gsub(indivT$taxon, pattern="p__", replacement="")

## Number of days and accum. degree days are strongly correlated.
with(indivT, cor(degdays, days))

## For each bacteria, plot counts for each pig vs. accum. degree days.
ggplot(indivT, aes(degdays, counts)) +
  geom_point() +
  facet_wrap(~species)

## Find species that have any day at which the count is over
## 479, which is the top 10% of counts.
commonSpecies <- indivT %>%
  filter(counts > 479) %>%
  select(species) %>%
  distinct(species) %>%
  unlist()

## For each bacteria, plot counts for each pig vs. accum. degree days.
ggplot(subset(indivT, species %in% commonSpecies), aes(degdays, counts)) +
  geom_point() +
  facet_wrap(~species)

## ##################################################





