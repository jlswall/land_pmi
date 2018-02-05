library("tidyverse")
library("randomForest")
library("figdim")

## Try using all the data, even though the observations at the end are
## taken less often than the ones at the beginning.


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0("../", taxalevel, "_massaged.csv"))
## ##################################################



## ##################################################
## Put the data in wide format.

## Move back to wide format.
wideT <- taxaT %>%
  ungroup() %>%
  select(days, degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay)
## ##################################################



## ##################################################
## Look at a cluster analysis, just to see whether observations made
## at the same time tend to be grouped together.  Results are mixed.
set.seed(4517802)
allT <- wideT %>% select(-subj, -Rare, -days)
hc.out <- hclust(dist(allT %>% select(-degdays)), method="average")
plot(hc.out, labels=allT$degdays)
## ##################################################



## ##################################################
## Try random forests for regression using "degdays" as the response
## variable.

## #########
## How many predictors?  (All columns except response: "degdays").
numPredictors <- ncol(allT) - 1

## Try different numbers of bootstrap samples.
numBtSampsVec <- seq(1000, 5000, by=1000)

## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
numVarSplitVec <- seq( floor(sqrt(numPredictors))-1, ceiling(numPredictors/3)+2 , by=2)
## numVarSplitVec <- seq( floor(sqrt(numPredictors)), ceiling(numPredictors/3)+3 , by=1)

## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)


## ## Uncomment the following to do cross-validation over and over.
## ## Number of times to do cross-validation.
## numCVs <- 100
## ## How many observations to reserve for testing each time.
## numLeaveOut <- round(0.10 * nrow(wideT))
## ## For matrix to hold cross-validation results.
## cvMSE <- matrix(NA, nrow(combos), ncol=numCVs)


## Number of times to do cross-validation.
numFolds <- 10
## Figure out which observations are in which fold by randomly
## assigning the numbers 1-10 to the various rows.  There are 93
## observations, so we assign an extra 1, 2, 3.
set.seed(305619)
numbersToAssign <- c( rep( 1:10, floor(nrow(allT)/10) ), 1:(nrow(allT) %% 10) )
whichFold <- sample(numbersToAssign, replace=F)

## For matrix to hold cross-validation results.
cvMSE <- matrix(NA, nrow(combos), ncol=numFolds)
cvErrFrac <- matrix(NA, nrow(combos), ncol=numFolds)
origUnitsqrtcvMSE <- matrix(NA, nrow(combos), ncol=numFolds)
origUnitsqrtcvErrFrac <- matrix(NA, nrow(combos), ncol=numFolds)
sqrtcvMSE <- matrix(NA, nrow(combos), ncol=numFolds)
sqrtcvErrFrac <- matrix(NA, nrow(combos), ncol=numFolds)


## Do cross-validation.
for (i in 1:numFolds){
  ## ############
  ## for (i in 1:numCVs){
  ## ##  Uncomment the following to do cross-validation over and over.
  ## ## Determine training and cross-validation set.
  ## whichLeaveOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)    
  ## subT <- wideT[-whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
  ## cvsetT <- wideT[whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
  ## ############

  ## Leave out all the rows assigned to i.
  whichLeaveOut <- whichFold==i
  cvsetT <- allT[whichLeaveOut,]
  subT <- allT[!whichLeaveOut,]

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
}
rm(i, j, SSTot)



combos$avgcvMSE <- apply(cvMSE, 1, mean)
combos$avgcvErrFrac <- apply(cvErrFrac, 1, mean)

combos$avgsqrtcvMSE <- apply(sqrtcvMSE, 1, mean)
combos$avgsqrtcvErrFrac <- apply(sqrtcvErrFrac, 1, mean)
combos$avgorigUnitsqrtcvMSE <- apply(origUnitsqrtcvMSE, 1, mean)
combos$avgorigUnitsqrtcvErrFrac <- apply(origUnitsqrtcvErrFrac, 1, mean)

write_csv(combos, path="all_data_avg_cv_metrics.csv")

## Exclude 14, 16, 18.  Try to figure out which lower number it should be.
ggplot(data=combos, aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line()
X11()
ggplot(data=combos, aes(x=numBtSamps, y=avgsqrtcvMSE, color=as.factor(numVarSplit))) + geom_line()
X11()
ggplot(data=combos, aes(x=numBtSamps, y=avgorigUnitsqrtcvMSE, color=as.factor(numVarSplit))) + geom_line()

## Eight, nine, or ten at each split looks best, with number of
## bootstrap samples indeterminate.
ggplot(data=combos, aes(x=numBtSamps, y=avgcvErrFrac, color=as.factor(numVarSplit))) + geom_line()
X11()
ggplot(data=combos, aes(x=numBtSamps, y=avgsqrtcvErrFrac, color=as.factor(numVarSplit))) + geom_line()
X11()
ggplot(data=combos, aes(x=numBtSamps, y=avgorigUnitsqrtcvErrFrac, color=as.factor(numVarSplit))) + geom_line()
## ##################################################
