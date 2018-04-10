library("tidyr")
library("readr")


## Read in data from Excel.
fileNm <- "raw_degdays_by_day_pig_from_baneshwar.csv"
rawT <- read_csv(file=fileNm)

## Remove the last line, which is a note giving the units as Celsius.
rawT <- rawT[-7,]
## Remove the last column, which is just the average degday days for
## each time step, across all pigs.
rawT <- rawT[,-7]


## Separate out and rename the first 2 columns for pig 1.  Add a column with the
## subject number.
pig1T <- rawT[,1:2]
colnames(pig1T) <- c("days", "degdays")
pig1T$subj <- "P1"
## Separate out and rename columns 3-4 for pig 2.  Add a column with the
## subject number.
pig2T <- rawT[,3:4]
colnames(pig2T) <- c("days", "degdays")
pig2T$subj <- "P2"
## Separate out and rename columns 5-6 for pig 3.  Add a column with the
## subject number.
pig3T <- rawT[,5:6]
colnames(pig3T) <- c("days", "degdays")
pig3T$subj <- "P3"


## Now, combined these into a long form dataset.
longT <- rbind(pig1T, pig2T, pig3T)
rm(pig1T, pig2T, pig3T)


## Remove the "T" from the "T0", "T1", etc. in the days column.
longT$days <- substring(longT$days, first=2, last=2)


write_csv(longT, path="degdays_by_subj_day.csv")
