library("readxl")
library("tidyr")


## library("openxlsx")
fileNm <- "../orig_data/Shane_all_skin_samples_taxo_bs_05_05_2017.xlsx"
## phylumAllT <- read.xlsx(fileNm, sheet="Phylum_all", skipEmptyCols=FALSE)
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
## of the observed counts over all pigs at the various times.

## We're ignoring columns past 114.
mainT <- rawAllT[,1:114]

## Identify column names starting with "A".
namesA <- colnames(mainT)[substring(first=1, last=1, colnames(mainT))=="A"]
wideIndivT <- mainT[,c("taxon", namesA)]
rm(namesA)

## Identify column names starting with "T".  These are the averages of
## the individual pigs at each time point.
namesT <- colnames(mainT)[substring(first=1, last=1, colnames(mainT))=="T"]
wideAvgsT <- mainT[,c("taxon", namesT)]
rm(namesT)

rm(mainT)
## ##################################################



## ##################################################
## Try to go from wide format to long format.

indivT <- wideIndivT %>%
  gather(indiv_time, counts, -taxon) %>%
  separate(indiv_time, into=c("subj", "time"))
## To see more info about how this works, see:
##   http://www.milanor.net/blog/reshape-data-r-tidyr-vs-reshape2/

## Next thing to check is whether we can get back the values in
## wideAvgsT when we do a summary for indivT.
## ##################################################

## Colums DN (phylum names, column 118) and DO-ED (columns 119-134)
## seem to be averages from the earlier columns.
phylumAllDF[1:5, 118:134]
## Column L matches column DO.
identical(phylumAllDF[,12], phylumAllDF[,119])
## Column S matches column DP.
identical(phylumAllDF[,19], phylumAllDF[,120])
## Column Z matches column DQ.
identical(phylumAllDF[,26], phylumAllDF[,121])
## Column DJ matches column ED.
identical(phylumAllDF[,114], phylumAllDF[,134])

## Names of phylums are almost same from column C to DN (118).  In
## column C, the names start with "p__", except for "unclassified" and
## "k__Bacteria".  So we remove "p__" and "k__" from those strings.
## We also capitalize "unclassifed" to that it will match
## "Unclassified) in column DN.
simplifiedColumnC <- gsub(pattern="p__", replacement="", phylumAllDF[,3] )
simplifiedColumnC <- gsub(pattern="k__", replacement="", simplifiedColumnC)
simplifiedColumnC <- gsub(pattern="unclassified", replacement="Unclassified", simplifiedColumnC)
identical(simplifiedColumnC, phylumAllDF[,118])

## Let's isolate these rows and columns.
smallDF <- phylumAllDF[,118:134]
## rownames(smallDF) <- phylumAllDF[,118]

library("ggplot2")
