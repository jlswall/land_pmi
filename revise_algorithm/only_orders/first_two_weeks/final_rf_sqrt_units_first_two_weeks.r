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

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0("../../", taxalevel, "_massaged.csv"), col_types="iiccnn")
## ##################################################



## ##################################################
## Put the data in wide format and restrict to the first 15 days.

## Move back to wide format.
earlyT <- taxaT %>%
  filter((days <= 15) & (taxa!="Rare")) %>%
  select(degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay) %>%
  select(-subj)

## Just for reference later, keep the days and degree days, so we can
## look at the time correspondence.
timeT <- taxaT %>% distinct(days, degdays)

rm(taxaT)
## ##################################################



## ##################################################
## From earlier experiments, we figured out these parameters work best
## for the random forest model.

## Number of bootstrap samples.
numBtSamps <- 5000

## Early runs indicated that the number of splits is around 16.
numVarSplit <- 16
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


set.seed(9402423)

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

set.seed(5850932)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(sqrt(degdays) ~ . , data=earlyT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

init.fig.dimen(file=paste0("sqrt_units_first_two_weeks_orders_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of order taxa (sqrt units, first 2 weeks)")
dev.off()


## In square root units:
## Find residuals:
resids <- rf$predicted - sqrt(earlyT$degdays)
## Print out RMSE:
sqrt( mean( resids^2 ) )
## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (sqrt(earlyT$degdays) - mean(sqrt(earlyT$degdays)))^2 ) )

## Projecting onto original units:
## Find estimated residuals:
resids <- (rf$predicted^2) - earlyT$degdays
## Print out RMSE:
sqrt( mean( resids^2 ) )
## Estimate of explained variance
1 - ( sum(resids^2)/sum( (earlyT$degdays - mean(earlyT$degdays))^2 ) )
## ##################################################
