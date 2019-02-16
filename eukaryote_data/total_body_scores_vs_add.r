library("tidyverse")
library("ggplot2")
library("readxl")

## Get the days and accumulated degree days from the other data we've
## already read in.
fileNm <- "families_massaged.csv"
rawdaysT <- read_csv(file=fileNm)
daysADDT <- unique(rawdaysT %>% select(days, degdays))
rm(rawdaysT, fileNm)

## Read in data from Excel.
fileNm <- "total_body_scores.xlsx"
listADDs <- as.vector(as.matrix(read_excel(path=fileNm, range="D1:S1", col_names=F)))

## Print out our mismatched ADDs and these, side by side.
write.table(file="mismatched_ADD_in_tbs_file.txt", cbind(c(daysADDT$degdays, NA), c(NA, listADDs)), sep="\t", row.names=F, col.names=F)


## ##################################################
## Code I was hoping to use to combine our ADDs with the total body
## scores doesn't work, because there are 19 days of total body scores
## and only 18 ADDs.
rawAllT <- read_excel(path=fileNm, range="A22:S27", col_names=F)
## Drop the column which contains just the string "TBS".
wideT <- rawAllT[,-2]
colnames(wideT) <- c("subj", paste("ADD.", daysADDT$degdays, sep=""))
rm(fileNm)

## Column names are unique in this sheet.
sum(duplicated(colnames(rawAllT)))

