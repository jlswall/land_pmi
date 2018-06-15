library("tidyverse")
library("randomForest")
library("figdim")
library("parallel")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0(taxalevel, "_only_first_two_weeks.csv"))
## ##################################################



## ##################################################
## Put the data in wide format; remove days, subj, and rare taxa.

## Move back to wide format.
wideT <- taxaT %>%
  filter(taxa!="Rare") %>%
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
numBtSamps <- 3000

## Using more cross-validations (1000) and leaving 20% out, it seems
## the number of variables to consider at each split should be about
## 12 (closely followed by 13).  This is using the response variable
## in sqrt units.
numVarSplit <- 12
## ##################################################



## ##################################################
## Run the cross-validation, for this model, so that we can see what
## the CV MSE looks like.

set.seed(7805903)

## Number of times to do cross-validation.
numCVs <- 1000
## How many observations to reserve for testing each time.
numLeaveOut <- round(0.20 * nrow(wideT))


## ###########################
## Set up function for fitting random forest model using square root
## units.
sqrtUnitsF <- function(x, mtry, ntree){
  sqrtrf <- randomForest(sqrt(degdays) ~ . , data=x$trainT, mtry=mtry, ntree=ntree, importance=T)
  return(predict(sqrtrf, newdata=x$validT))
}
## ###########################


## ###########################
## Get set up for cross-validation.
crossvalidL <- vector("list", numCVs)
for (i in 1:numCVs){
  lvOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)
  trainT <- wideT[-lvOut,]
  validT <- wideT[lvOut,]
  crossvalidL[[i]] <- list(trainT=trainT, validT=validT)
}
rm(i, lvOut, trainT, validT)

## Conduct cross-validation.
sqrtFitL <- mclapply(crossvalidL, mc.cores=6, sqrtUnitsF, mtry=numVarSplit, ntree=numBtSamps)
## ###########################


## ###########################
## For matrix to hold cross-validation results.
sqrtcvMSE <- rep(NA, numCVs)
sqrtcvErrFrac <- rep(NA, numCVs)
origUnitsqrtcvMSE <- rep(NA, numCVs)
origUnitsqrtcvErrFrac <- rep(NA, numCVs)

set.seed(5432472)

## Now, calculate the various summary statistics for each cross-validation.
for (i in 1:numCVs){

  ## Get the validation set for this run from the list.
  validT <- crossvalidL[[i]][["validT"]]

  ## Calculate SSTotal for the cross-validation set.
  SSTot <- sum( (validT$degdays-mean(validT$degdays))^2 )

  ## Calculate the MSE and error fraction of the SS Total for the
  ## validation data in the original units.
  sqrtUnitResid <- sqrtFitL[[i]] - sqrt(validT$degdays)
  origUnitResid <- sqrtFitL[[i]]^2 - validT$degdays
  sqrtcvMSE[i] <- mean(sqrtUnitResid^2)
  sqrtcvErrFrac[i] <- sum(sqrtUnitResid^2) / sum( ( sqrt(validT$degdays) - mean(sqrt(validT$degdays)) )^2 )
  origUnitsqrtcvMSE[i] <- mean(origUnitResid^2)
  origUnitsqrtcvErrFrac[i] <- sum(origUnitResid^2)/SSTot
  rm(sqrtUnitResid, origUnitResid)
}
rm(i, validT, SSTot)

write_csv(data.frame(sqrtcvMSE, sqrtcvErrFrac, origUnitsqrtcvMSE, origUnitsqrtcvErrFrac), path="final_rf_sqrt_units_cvstats_first_two_weeks.csv")
rm(sqrtcvMSE, sqrtcvErrFrac, origUnitsqrtcvMSE, origUnitsqrtcvErrFrac)
## ##################################################



## ##################################################
## Fit the final random forest with all the data (no cross-validation).

set.seed(3817139)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(sqrt(degdays) ~ . , data=wideT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

init.fig.dimen(file=paste0("sqrt_units_first_two_weeks_families_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of family taxa (sqrt. units, first two weeks)")
dev.off()


## In square root units:
## Find residuals:
resids <- rf$predicted - sqrt(wideT$degdays)

## Print out RMSE:
sqrt( mean( resids^2 ) )
## RMSE: 2.307899

## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (sqrt(wideT$degdays) - mean(sqrt(wideT$degdays)))^2 ) )
## Expl. frac.: 0.8786743


## Projecting onto original units:
## Find estimated residuals:
resids <- (rf$predicted^2) - wideT$degdays

## Print out RMSE:
sqrt( mean( resids^2 ) )
## RMSE projected to original units: 64.09809

## Estimate of explained variance
1 - ( sum(resids^2)/sum( (wideT$degdays - mean(wideT$degdays))^2 ) )
## Expl. frac. projected to orig. units: 0.8247467
## ##################################################
