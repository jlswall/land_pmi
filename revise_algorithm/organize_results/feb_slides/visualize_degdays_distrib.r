library("tidyverse")


## Read in the orders data just to get the distribution of degree
## days.
taxaT <- read_csv("../families_massaged.csv", col_types="iiccnn")
## Just for reference later, keep the days and degree days, so we can
## look at the time correspondence.
timeT <- taxaT %>% distinct(days, degdays)
rm(taxaT)


## ##########
## Make histograms of degree days (orig. and sqrt. units).

## Histogram of degree days (all time steps).
pdf(file="degdays_all_time_steps_hist.pdf")
hist(timeT$degdays, xlab="Degree days", main=NULL)
dev.off()

## Histogram of degree days (all time steps).
pdf(file="sqrt_degdays_all_time_steps_hist.pdf")
hist(sqrt(timeT$degdays), xlab="Square root of degree days", main=NULL)
dev.off()
## ##########


## ##########
## Plot only the first 15 days.

subT <- timeT %>% filter(days <= 15)

## Histogram of degree days (first 15 days)
pdf(file="degdays_first_two_weeks_hist.pdf")
hist(subT$degdays, xlab="Degree days (first 15 days)", main=NULL)
dev.off()

## Histogram of degree days (first 15 days)
pdf(file="sqrt_degdays_first_two_weeks_hist.pdf")
hist(sqrt(subT$degdays), xlab="Square root of degree days (first 15 days)", main=NULL)
dev.off()
## ##########
