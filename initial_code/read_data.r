library("openxlsx")
fileNm <- "../orig_data/Shane_all_skin_samples_taxo_bs_05_05_2017.xlsx"
phylumAllDF <- read.xlsx(fileNm, sheet="Phylum_all", skipEmptyCols=FALSE)

## Some of these column names are duplicated.
duplicated(colnames(phylumAllDF))
## This first happens at column 116.
min(which(duplicated(colnames(phylumAllDF))))


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
