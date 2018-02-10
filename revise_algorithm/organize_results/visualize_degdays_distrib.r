library("tidyverse")
library("figdim")


## Read in the orders data just to get the distribution of degree
## days.
taxaT <- read_csv("../families_massaged.csv", col_types="iiccnn")
## Just for reference later, keep the days and degree days, so we can
## look at the time correspondence.
timeT <- taxaT %>% distinct(days, degdays)
rm(taxaT)


## Histogram of degree days.
pdf(file="hist_degdays.pdf")
hist(timeT$degdays, xlab="Degree days", main=NULL)
dev.off()

## Histogram of degree days.
pdf(file="hist_sqrt_degdays.pdf")
hist(sqrt(timeT$degdays), xlab="Square root of degree days", main=NULL)
dev.off()


