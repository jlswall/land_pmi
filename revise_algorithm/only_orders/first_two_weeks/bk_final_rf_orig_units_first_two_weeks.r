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
taxaT <- read_csv(paste0("../../", taxalevel, "_massaged.csv"))
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

## Early runs indicated that the number of variables to consider at
## each split is about 8 for response variable in the original units.
numVarSplit <- 8
## ##################################################


## ##################################################
## Run the cross-validation for this model, so that we can see what
## the CV MSE looks like.

## Number of times to do cross-validation.
numCVs <- 100
## ## How many observations to reserve for testing each time.
numLeaveOut <- round(0.10 * nrow(earlyT))


## For matrix to hold cross-validation results.
cvMSE <- rep(NA, numCVs)
cvErrFrac <- rep(NA, numCVs)

set.seed(943179)


## Do cross-validation.
for (i in 1:numCVs){
  
  ## Determine training and cross-validation set.
  whichLeaveOut <- sample(1:nrow(earlyT), size=numLeaveOut, replace=F)    
  subT <- earlyT[-whichLeaveOut,]
  cvsetT <- earlyT[whichLeaveOut,]
  
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

write_csv(data.frame(cvMSE, cvErrFrac), path="final_rf_orig_units_cvstats_first_two_weeks.csv")
rm(cvMSE, cvErrFrac)
## ##################################################



## ##################################################
## Fit the final random forest with all the data (no cross-validation).

set.seed(880932)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(degdays ~ . , data=earlyT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

init.fig.dimen(file=paste0("orig_units_first_two_weeks_orders_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of order taxa (orig. units, first 2 weeks)")
dev.off()


## Find residuals:
resids <- rf$predicted - earlyT$degdays
## Print out RMSE:
sqrt( mean( resids^2 ) )
## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (earlyT$degdays - mean(earlyT$degdays))^2 ) )
## ##################################################
