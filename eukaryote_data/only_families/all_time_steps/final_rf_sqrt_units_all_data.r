library("tidyverse")
library("randomForest")
library("figdim")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0("../../", taxalevel, "_massaged.csv"))
## ##################################################



## ##################################################
## Put the data in wide format; remove days, subj, and rare taxa.

## Move back to wide format.
allT <- taxaT %>%
  filter(taxon!="Rare") %>%
  select(degdays, subj, taxon, fracBySubjDay) %>%
  spread(taxon, fracBySubjDay) %>%
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
numBtSamps <- 3000

## Using more cross-validations and leaving 20% out, the number of
## variables to consider at each split should be around 17 (with 16
## and 18 very close), modeling the response variable in sqrt units.
numVarSplit <- 17
## ##################################################



## ##################################################
## Run the cross-validation for this model, so that we can see what
## the CV MSE looks like.

## Number of times to do cross-validation.
numCVs <- 500
## ## How many observations to reserve for testing each time.
numLeaveOut <- round(0.20 * nrow(allT))

## For matrix to hold cross-validation results.
sqrtcvMSE <- rep(NA, numCVs)
sqrtcvErrFrac <- rep(NA, numCVs)
origUnitsqrtcvMSE <- rep(NA, numCVs)
origUnitsqrtcvErrFrac <- rep(NA, numCVs)


set.seed(3450460)

## Do cross-validation.
for (i in 1:numCVs){
  
  ## Determine training and cross-validation set.
  whichLeaveOut <- sample(1:nrow(allT), size=numLeaveOut, replace=F)    
  subT <- allT[-whichLeaveOut,]
  cvsetT <- allT[whichLeaveOut,]
  
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


write_csv(data.frame(sqrtcvMSE, sqrtcvErrFrac, origUnitsqrtcvMSE, origUnitsqrtcvErrFrac), path="final_rf_sqrt_units_cvstats_all_data.csv")
rm(sqrtcvMSE, sqrtcvErrFrac, origUnitsqrtcvMSE, origUnitsqrtcvErrFrac)
## ##################################################



## ##################################################
## Fit the final random forest with all the data (no cross-validation).

set.seed(5057430)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(sqrt(degdays) ~ . , data=allT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

init.fig.dimen(file=paste0("sqrt_units_all_data_famliess_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of family taxa (sqrt. units, all time steps)")
dev.off()


## In square root units:
## Find residuals:
resids <- rf$predicted - sqrt(allT$degdays)
## Print out RMSE:
sqrt( mean( resids^2 ) )
## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (sqrt(allT$degdays) - mean(sqrt(allT$degdays)))^2 ) )

## Projecting onto original units:
## Find estimated residuals:
resids <- (rf$predicted^2) - allT$degdays
## Print out RMSE:
sqrt( mean( resids^2 ) )
## Estimate of explained variance
1 - ( sum(resids^2)/sum( (allT$degdays - mean(allT$degdays))^2 ) )
## ##################################################
