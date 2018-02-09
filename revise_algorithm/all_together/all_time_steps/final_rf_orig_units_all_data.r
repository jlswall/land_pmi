library("tidyverse")
library("randomForest")
library("figdim")


## ##################################################
## Read in the data.

## Read in the combined data in wide format.
tmpT <- read_csv("combined_taxa.csv")

## Remove subj and days variable, since they won't be used in the
## random forest.
allT <- tmpT %>%
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

## Early runs indicated that the number of splits is around 45-50.
numVarSplit <- 45
## ##################################################



## ##################################################
## Run the cross-validation for this model, so that we can see what
## the CV MSE looks like.

## Number of times to do cross-validation.
numCVs <- 100
## ## How many observations to reserve for testing each time.
numLeaveOut <- round(0.10 * nrow(allT))

## For matrix to hold cross-validation results.
cvMSE <- rep(NA, numCVs)
cvErrFrac <- rep(NA, numCVs)

## Do cross-validation.
for (i in 1:numCVs){
  
  ## Determine training and cross-validation set.
  whichLeaveOut <- sample(1:nrow(allT), size=numLeaveOut, replace=F)    
  subT <- allT[-whichLeaveOut,]
  cvsetT <- allT[whichLeaveOut,]
  
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


write_csv(data.frame(cvMSE, cvErrFrac), path="final_rf_cvstats_all_data.csv")
rm(cvMSE, cvErrFrac)
## ##################################################



## ##################################################
## Fit the final random forest with all the data (no cross-validation).

set.seed(603932)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(degdays ~ . , data=allT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

init.fig.dimen(file=paste0("orig_units_combined_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of combined taxa (all time steps)")
dev.off()
