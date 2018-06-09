library("tidyverse")
library("randomForest")
library("figdim")
library("parallel")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "orders"

## Read in cleaned-up phyla, orders, or families taxa for JUST THE
## FIRST 2 WEEKS.
earlyT <- read_csv(paste0(taxalevel, "_only_first_two_weeks.csv"))
## ##################################################




## ##################################################
## Put the data in wide format; remove days, subj, and rare taxa.

## Move back to wide format.
wideT <- earlyT %>%
  filter(taxon!="Rare") %>%
  select(degdays, subj, taxon, fracBySubjDay) %>%
  spread(taxon, fracBySubjDay) %>%
  select(-subj)
## ##################################################



## ##################################################
## Try random forests for regression using "degdays" as the response
## variable.

## #########
## How many predictors?  (All columns except response: "degdays").
numPredictors <- ncol(wideT) - 1

## Try different numbers of bootstrap samples.
numBtSampsVec <- c(900, 1500, 2100, 2400, 3000)

## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
numVarSplitVec <- seq(15, 33, by=1)

## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)


## ###########################
## Do cross-validation over and over, leaving out a different 20% of
## the observations each time.

set.seed(9838462)

## Number of times to do cross-validation.
numCVs <- 1000
## How many observations to reserve for testing each time.
numLeaveOut <- round(0.20 * nrow(wideT))


## For matrix to hold cross-validation results.
cvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
cvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)
origUnitsqrtcvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
origUnitsqrtcvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)
sqrtcvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
sqrtcvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)



## #########################################
## Set up function for fitting random forest model using original
## units.
origUnitsF <- function(x, jCombo){
  rf <- randomForest(degdays ~ . , data=x$trainT, mtry=combos[jCombo, "numVarSplit"], ntree=combos[jCombo, "numBtSamps"], importance=T)
  return(predict(rf, newdata=x$validT))
}

## Set up function for fitting random forest model using square root
## units.
sqrtUnitsF <- function(x, jCombo){
  sqrtrf <- randomForest(sqrt(degdays) ~ . , data=x$trainT, mtry=combos[jCombo, "numVarSplit"], ntree=combos[jCombo, "numBtSamps"], importance=T)
  return(predict(sqrtrf, newdata=x$validT))
}
## #########################################


## #########################################
## For each cross-validation we plan to do, randomly select training
## and validation training sets.
crossvalidL <- vector("list", numCVs)
for (i in 1:numCVs){
  lvOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)
  trainT <- wideT[-lvOut,]
  validT <- wideT[lvOut,]
  crossvalidL[[i]] <- list(trainT=trainT, validT=validT)
}
rm(i, lvOut, trainT, validT)
## #########################################


## #########################################
## Try using lapply to fit the random forests.
origFitL <- vector("list", nrow(combos))
for (j in 1:nrow(combos)){
  origFitL[[j]] <- mclapply(crossvalidL, mc.cores=6, origUnitsF, jCombo=j)
  if (j %% 2 == 0)
    print(paste0("In orig units, finished combo number ", j))
}    
rm(j)

sqrtFitL <- vector("list", nrow(combos))
for (j in 1:nrow(combos)){
  sqrtFitL[[j]] <- mclapply(crossvalidL, mc.cores=6, sqrtUnitsF, jCombo=j)
  if (j %% 2 == 0)
    print(paste0("In sqrt units, finished combo number ", j))
}
rm(j)
## #########################################


## #########################################
## Now, calculate the various summary statistics for each cross-validation.

for (i in 1:numCVs){

  ## Get the validation set for this run from the list.
  validT <- crossvalidL[[i]][["validT"]]

  ## Calculate SSTotal for the cross-validation set.
  SSTot <- sum( (validT$degdays-mean(validT$degdays))^2 )
  
  for (j in 1:nrow(combos)){
     
    ## Calculate the MSE and error fraction of the SS Total for the
    ## validation data in the original units.
    resid <- origFitL[[j]][[i]] - validT$degdays
    cvMSE[j,i] <- mean(resid^2)
    cvErrFrac[j,i] <- sum(resid^2)/SSTot
    rm(resid)
  
    ## Calculate the MSE and error fraction of the SS Total for the
    ## validation data in the original units.
    sqrtUnitResid <- sqrtFitL[[j]][[i]] - sqrt(validT$degdays)
    origUnitResid <- sqrtFitL[[j]][[i]]^2 - validT$degdays
    sqrtcvMSE[j,i] <- mean(sqrtUnitResid^2)
    sqrtcvErrFrac[j,i] <- sum(sqrtUnitResid^2)/sum( ( sqrt(validT$degdays) - mean(sqrt(validT$degdays)) )^2 )
    origUnitsqrtcvMSE[j,i] <- mean(origUnitResid^2)
    origUnitsqrtcvErrFrac[j,i] <- sum(origUnitResid^2)/SSTot
    rm(sqrtUnitResid, origUnitResid)
  }
}
rm(i, j, validT, SSTot)
## #########################################




combos$avgcvMSE <- apply(cvMSE, 1, mean)
combos$avgcvErrFrac <- apply(cvErrFrac, 1, mean)

combos$avgsqrtcvMSE <- apply(sqrtcvMSE, 1, mean)
combos$avgsqrtcvErrFrac <- apply(sqrtcvErrFrac, 1, mean)
combos$avgorigUnitsqrtcvMSE <- apply(origUnitsqrtcvMSE, 1, mean)
combos$avgorigUnitsqrtcvErrFrac <- apply(origUnitsqrtcvErrFrac, 1, mean)

write_csv(combos, path="parallel_leave_out_20perc_first_two_weeks.csv")


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
