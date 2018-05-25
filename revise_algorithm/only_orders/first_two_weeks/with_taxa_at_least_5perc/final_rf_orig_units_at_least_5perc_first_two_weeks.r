library("tidyverse")
library("randomForest")
library("figdim")

## In the first 15 days (448 degree days), observations were made
## approx. every other day.  Then, there's a big gap between 15 days
## and the next observation 26 days out.  Let's try the analysis with
## just these observations.


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "orders"

## Read in cleaned-up phyla, orders, or families taxa for JUST THE
## FIRST 2 WEEKS, using only taxa with at least %% on at least one
## subject and day during this time period.
earlyT <- read_csv(paste0(taxalevel, "_at_least_5perc_only_first_two_weeks.csv"))
## ##################################################



## ##################################################
## Put the data in wide format and restrict to the first 15 days.

## Move to wide format.
wideT <- earlyT %>%
  filter(taxa!="Rare") %>%
  select(degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay) %>%
  select(-subj)

## Just for reference later, keep the days and degree days, so we can
## look at the time correspondence.
timeT <- earlyT %>% distinct(days, degdays)

rm(earlyT)
## ##################################################



## ##################################################
## From earlier experiments, we figured out these parameters work best
## for the random forest model.

## Number of bootstrap samples.
numBtSamps <- 5000

## Simulation runs (with 20% of data reserved for validation)
## indicated that the number of variables to consider at each split is
## 7 (closely followed by 8) for response variable in the original
## units.  Interestingly, this is the same number chosen when we
## allowed taxa making up at least 1%.
numVarSplit <- 7
## ##################################################


## ##################################################
## Run the cross-validation for this model, so that we can see what
## the CV MSE looks like.

## Number of times to do cross-validation.
numCVs <- 100
## How many observations to reserve for testing each time.
numLeaveOut <- round(0.20 * nrow(wideT))


## For matrix to hold cross-validation results.
cvMSE <- rep(NA, numCVs)
cvErrFrac <- rep(NA, numCVs)

set.seed(923407)


## Do cross-validation.
for (i in 1:numCVs){
  
  ## Determine training and cross-validation set.
  whichLeaveOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)    
  subT <- wideT[-whichLeaveOut,]
  cvsetT <- wideT[whichLeaveOut,]
  
  ## Calculate SSTotal for the cross-validation set.
  SSTot <- sum( (cvsetT$degdays-mean(cvsetT$degdays))^2 )

  rf <- randomForest(degdays ~ . , data=subT, mtry=numVarSplit,
                     ntree=numBtSamps, importance=T)
  fitTest <- predict(rf, newdata=cvsetT)
    
  fitResid <- fitTest - cvsetT$degdays
  cvMSE[i] <- mean(fitResid^2)
  cvErrFrac[i] <- sum(fitResid^2)/SSTot
}
rm(whichLeaveOut, subT, cvsetT, SSTot, rf, fitTest, fitResid, numLeaveOut, i)

write_csv(data.frame(cvMSE, cvErrFrac), path="final_rf_orig_units_cvstats_at_least_5perc_first_two_weeks.csv")
rm(cvMSE, cvErrFrac)
## ##################################################



## ##################################################
## Fit the final random forest with all the data (no cross-validation).

set.seed(5384432)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(degdays ~ . , data=wideT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

init.fig.dimen(file=paste0("orig_units_at_least_5perc_first_two_weeks_orders_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of order taxa (orig. units, first 2 weeks)")
dev.off()


## Find residuals:
resids <- rf$predicted - wideT$degdays
## Print out RMSE:
sqrt( mean( resids^2 ) )
## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (wideT$degdays - mean(wideT$degdays))^2 ) )
## ##################################################
