library("tidyverse")
library("randomForest")
library("figdim")

## In the first 15 days (448 degree days), observations were made
## approx. every other day.  Then, there's a big gap between 15 days
## and the next observation 26 days out.  Let's try the analysis with
## just these observations.


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "phyla"

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
## Try random forests for regression using "days" as the response
## variable.

## #########
## How many predictors?  (All columns except response: "degdays").
numPredictors <- ncol(earlyT) - 1

## Try different numbers of bootstrap samples.
numBtSampsVec <- seq(4000, 5000, by=1000)

## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
## Our first run showed that 5 was best for the square root response,
## but on theoriginal units, 6-10 were grouped tightly together.  So,
## I'm doing another run with more cross-validation samples, different
## random seed, and excluding 2 and 3.
numVarSplitVec <- seq(4, numPredictors, by=1)

## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)


## ###########################
## Do cross-validation over and over, leaving out a different 10% of
## the 57 observations each time.

set.seed(9989932)

## Number of times to do cross-validation.
numCVs <- 200
## ## How many observations to reserve for testing each time.
numLeaveOut <- round(0.10 * nrow(earlyT))


## For matrix to hold cross-validation results.
cvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
cvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)
origUnitsqrtcvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
origUnitsqrtcvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)
sqrtcvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
sqrtcvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)


## Do cross-validation.
for (i in 1:numCVs){
  
  ## Determine training and cross-validation set.
  whichLeaveOut <- sample(1:nrow(earlyT), size=numLeaveOut, replace=F)    
  subT <- earlyT[-whichLeaveOut,]
  cvsetT <- earlyT[whichLeaveOut,]
  
  ## Calculate SSTotal for the cross-validation set.
  SSTot <- sum( (cvsetT$degdays-mean(cvsetT$degdays))^2 )

  for (j in 1:nrow(combos)){    
    rf <- randomForest(degdays ~ . , data=subT, mtry=combos[j,"numVarSplit"], ntree=combos[j,"numBtSamps"], importance=T)
    fitTest <- predict(rf, newdata=cvsetT)
    
    fitResid <- fitTest - cvsetT$degdays
    cvMSE[j,i] <- mean(fitResid^2)
    cvErrFrac[j,i] <- sum(fitResid^2)/SSTot
  }
  rm(rf, fitTest, fitResid)

  ## In sqrt units:
  for (j in 1:nrow(combos)){    
    sqrtrf <- randomForest(sqrt(degdays) ~ . , data=subT, mtry=combos[j,"numVarSplit"], ntree=combos[j,"numBtSamps"], importance=T)
    sqrtfitTest <- predict(sqrtrf, newdata=cvsetT)
    sqrtfitResid <- sqrtfitTest - sqrt(cvsetT$degdays)
    origUnitResid <- sqrtfitTest^2 - cvsetT$degdays
    
    sqrtcvMSE[j,i] <- mean(sqrtfitResid^2)
    sqrtcvErrFrac[j,i] <- sum(sqrtfitResid^2)/sum( ( sqrt(cvsetT$degdays) - mean(sqrt(cvsetT$degdays)) )^2 )
    origUnitsqrtcvMSE[j,i] <- mean(origUnitResid^2)
    origUnitsqrtcvErrFrac[j,i] <- sum(origUnitResid^2)/SSTot
  }
  rm(sqrtrf, sqrtfitTest, sqrtfitResid, origUnitResid)

  if (i %% 10 == 0)
    print(paste0("Finishing cross-validation number ", i))
}
rm(i, j, SSTot)



combos$avgcvMSE <- apply(cvMSE, 1, mean)
combos$avgcvErrFrac <- apply(cvErrFrac, 1, mean)

combos$avgsqrtcvMSE <- apply(sqrtcvMSE, 1, mean)
combos$avgsqrtcvErrFrac <- apply(sqrtcvErrFrac, 1, mean)
combos$avgorigUnitsqrtcvMSE <- apply(origUnitsqrtcvMSE, 1, mean)
combos$avgorigUnitsqrtcvErrFrac <- apply(origUnitsqrtcvErrFrac, 1, mean)

write_csv(combos, path="repeated_cv_first_two_weeks_avg_cv_metrics.csv")


ggplot(data=combos, aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line()
## X11()
ggplot(data=combos, aes(x=numBtSamps, y=avgsqrtcvMSE, color=as.factor(numVarSplit))) + geom_line()
## X11()
ggplot(data=combos, aes(x=numBtSamps, y=avgorigUnitsqrtcvMSE, color=as.factor(numVarSplit))) + geom_line()


ggplot(data=combos, aes(x=numBtSamps, y=avgcvErrFrac, color=as.factor(numVarSplit))) + geom_line()
## X11()
ggplot(data=combos, aes(x=numBtSamps, y=avgsqrtcvErrFrac, color=as.factor(numVarSplit))) + geom_line()
## X11()
ggplot(data=combos, aes(x=numBtSamps, y=avgorigUnitsqrtcvErrFrac, color=as.factor(numVarSplit))) + geom_line()
## ####################
