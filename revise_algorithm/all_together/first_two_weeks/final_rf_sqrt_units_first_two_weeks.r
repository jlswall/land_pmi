library("tidyverse")
library("randomForest")
library("figdim")


## In the first 15 days (448 degree days), observations were made
## approx. every other day.  Then, there's a big gap between 15 days
## and the next observation 26 days out.  Let's try the analysis with
## just these observations.


## ##################################################
## Read in the data.

## Read in the combined data in wide format.
tmpT <- read_csv("../combined_taxa.csv")

## Remove subj and days variable, since they won't be used in the
## random forest.
earlyT <- tmpT %>%
  filter(days <= 15) %>%
  select(-subj, -days)

## Just for reference later, keep the days and degree days, so we can
## look at the time correspondence.
timeT <- tmpT %>% distinct(days, degdays)

rm(tmpT)
## ##################################################



## ##################################################
## From earlier experiments, we figured out these parameters work best
## for the random forest model.

## Number of bootstrap samples.
numBtSamps <- 5000

## Early runs indicated that the number of splits is around 25.
numVarSplit <- 25
## ##################################################



## ##################################################
## Run the cross-validation for this model, so that we can see what
## the CV MSE looks like.

## Number of times to do cross-validation.
numCVs <- 100
## ## How many observations to reserve for testing each time.
numLeaveOut <- round(0.10 * nrow(earlyT))

## For matrix to hold cross-validation results.
sqrtcvMSE <- rep(NA, numCVs)
sqrtcvErrFrac <- rep(NA, numCVs)
origUnitsqrtcvMSE <- rep(NA, numCVs)
origUnitsqrtcvErrFrac <- rep(NA, numCVs)

set.seed(4012423)

## Do cross-validation.
for (i in 1:numCVs){
  
  ## Determine training and cross-validation set.
  whichLeaveOut <- sample(1:nrow(earlyT), size=numLeaveOut, replace=F)    
  subT <- earlyT[-whichLeaveOut,]
  cvsetT <- earlyT[whichLeaveOut,]
  
  ## Calculate SSTotal for the cross-validation set.
  SSTot <- sum( (cvsetT$degdays-mean(cvsetT$degdays))^2 )

  sqrtrf <- randomForest(sqrt(degdays) ~ . , data=subT, mtry=numVarSplit, ntree=numBtSamps, importance=T)
  sqrtfitTest <- predict(sqrtrf, newdata=cvsetT)
  sqrtfitResid <- sqrtfitTest - sqrt(cvsetT$degdays)
  origUnitResid <- sqrtfitTest^2 - cvsetT$degdays
    
  sqrtcvMSE[i] <- mean(sqrtfitResid^2)
  sqrtcvErrFrac[i] <- sum(sqrtfitResid^2)/sum( ( sqrt(cvsetT$degdays) - mean(sqrt(cvsetT$degdays)) )^2 )
  origUnitsqrtcvMSE[i] <- mean(origUnitResid^2)
  origUnitsqrtcvErrFrac[i] <- sum(origUnitResid^2)/SSTot
}
rm(whichLeaveOut, subT, cvsetT, SSTot, i)
rm(sqrtrf, sqrtfitTest, sqrtfitResid, origUnitResid)


write_csv(data.frame(sqrtcvMSE, sqrtcvErrFrac, origUnitsqrtcvMSE, origUnitsqrtcvErrFrac), path="final_rf_sqrt_units_cvstats_first_two_weeks.csv")
rm(sqrtcvMSE, sqrtcvErrFrac, origUnitsqrtcvMSE, origUnitsqrtcvErrFrac)
## ##################################################



## ##################################################
## Fit the final random forest with the first two weeks of the data
## (no cross-validation).

set.seed(480932)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(sqrt(degdays) ~ . , data=earlyT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

init.fig.dimen(file=paste0("sqrt_units_first_two_weeks_combined_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of combined taxa (sqrt units, first 2 weeks)")
dev.off()
## ##################################################
