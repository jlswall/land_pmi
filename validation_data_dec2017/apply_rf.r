library("tidyverse")
library("randomForest")
library("readr")

## Load the random forest model that we fit using the first 2 weeks of
## the original data.
load("../revise_algorithm/only_orders/first_two_weeks/final_rf_orders_first_two_weeks.RData")

## Get list of variables considered in the model.
rfCovar <- attr(rf[["terms"]], "term.labels")


## Read in the validation data.
validT <- read_csv("orders_massaged.csv")


## For the validation data, assess taxa variability across the first
## few days, averaged across pigs, in terms of fractions.
avgTaxaOverdegdaysT <- validT %>%
  filter(degdays <= 150) %>%
  group_by(degdays, taxa) %>%
  summarize(avgFracByDegday=mean(fracBySubjDay)) %>%
  filter(avgFracByDegday >= 0.01)


ggplot(validT %>% filter(taxa %in% rfCovar),
    aes(degdays, fracBySubjDay)) +
    geom_point(aes(color=taxa)) +
    labs(x="Degree days", y="Composition fraction") +
    labs(color="Taxa")

ggplot(avgTaxaOverdegdaysT, aes(x=degdays, y=avgFracByDegday, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  labs(x="Days", y="Composition fraction")
         



## ##################################################
## Put the data in wide format.

## Move back to wide format.
wideT <- validT %>%
  select(degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay) %>%
  select(-subj)


## There's one taxa from the original model which isn't in the
## validation data.  That is "Streptophyta".
rfCovar[!(rfCovar %in% colnames(wideT))]
## Add a column of this name with zeroes in it.
wideT$Streptophyta <- 0.0
## ##################################################


## Make prediction for the validation set, using the original model.
predictValid <- predict(rf, newdata=wideT)
## Find residuals:
resids <- predictValid - wideT$degdays
## Print out RMSE:
sqrt( mean( resids^2 ) )
## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (wideT$degdays - mean(wideT$degdays))^2 ) )

## The predictions are not good, and I think this is because the
## fractions for the various taxa are quite different.


