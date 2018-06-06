library("tidyverse")
library("randomForest")
library("figdim")
library("parallel")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "orders"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0(taxalevel, "_at_least_5perc_all_time_steps.csv"))
## ##################################################



## ##################################################
## Put the data in wide format; remove days and subj.
## This time, LEAVE IN rare taxa.

## Move back to wide format.
allT <- taxaT %>%
  select(degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay) %>%
  select(-subj)

## Just for reference later, keep the days and degree days, so we can
## look at the time correspondence.
timeT <- taxaT %>% distinct(days, degdays)

rm(taxaT)
## ##################################################



## ##################################################
## Set up the parameters that govern our cross-validation runs for
## random forest models.

## How many predictors?  (All columns except response: "degdays").
numPredictors <- ncol(allT) - 1

## Try different numbers of bootstrap samples.
numBtSampsVec <- c(600, 1800, 2400, 3000)

## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
numVarSplitVec <- seq(6, 10, by=1)

## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)

rm(numBtSampsVec, numVarSplitVec)
## ##################################################



## ##################################################
## Set up functions for fitting random forest model using original
## and sqrt units.
origUnitsF <- function(x, jCombo){
  rf <- randomForest(degdays ~ . , data=x$trainT, mtry=combos[jCombo, "numVarSplit"], ntree=combos[jCombo, "numBtSamps"], importance=T)
  return(predict(rf, newdata=x$validT))
}

## ## Set up function for fitting random forest model using square root
## ## units.
## sqrtUnitsF <- function(x, jCombo){
##   sqrtrf <- randomForest(sqrt(degdays) ~ . , data=x$trainT, mtry=combos[jCombo, "numVarSplit"], ntree=combos[jCombo, "numBtSamps"], importance=T)
##   return(predict(sqrtrf, newdata=x$validT))
## }
## ##################################################




## ##################################################
## CROSS-VALIDATION RUNS: Randomly choose 20% of degdays to leave out
## in each run.  This means certain degdays will be completely missing
## in the training dataset.
## There are 16 time steps, so 20% means leave out about 3 days each
## time.  This means there are (16 choose 3) = 560 different
## cross-validation combinations.


## ####################
## Form matrix in which each row gives the 3 degree days to leave out.
allSubsetCombos <- t(combn(timeT$degdays, m=3))
numCVs <- nrow(allSubsetCombos)
## Set up training and validation sets using the sequence of
## combinations to leave out.
crossvalidL <- vector("list", numCVs)
for (i in 1:numCVs){
  trainT <- allT %>% filter(!(degdays %in% allSubsetCombos[i,]))
  validT <- allT %>% filter(degdays %in% allSubsetCombos[i,])
  crossvalidL[[i]] <- list(trainT=trainT, validT=validT)
}
rm(i, trainT, validT)
## ####################


## ####################
## Try using lapply to fit the random forests.

set.seed(6856965)

origFitL <- vector("list", nrow(combos))
for (j in 1:nrow(combos)){
  origFitL[[j]] <- mclapply(crossvalidL, mc.cores=7, origUnitsF, jCombo=j)
  if (j %% 2 == 0)
    print(paste0("Leaving out groups of days, finished combo number ", j))
}    
rm(j)
## ####################


## #####################
## Form matrices to hold cross-validation results.
cvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
cvErrFrac <- matrix(NA, nrow(combos), ncol=numCVs)
## #####################


## #####################
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
  
    ## ## Calculate the MSE and error fraction of the SS Total for the
    ## ## validation data in the original units.
    ## sqrtUnitResid <- sqrtFitL[[j]][[i]] - sqrt(validT$degdays)
    ## origUnitResid <- sqrtFitL[[j]][[i]]^2 - validT$degdays
    ## sqrtcvMSE[j,i] <- mean(sqrtUnitResid^2)
    ## sqrtcvErrFrac[j,i] <- sum(sqrtUnitResid^2)/sum( ( sqrt(validT$degdays) - mean(sqrt(validT$degdays)) )^2 )
    ## origUnitsqrtcvMSE[j,i] <- mean(origUnitResid^2)
    ## origUnitsqrtcvErrFrac[j,i] <- sum(origUnitResid^2)/SSTot
    ## rm(sqrtUnitResid, origUnitResid)
  }
}
rm(i, j, validT, SSTot)

## Average these summary statistics over all cross-validation runs.
## Original units:
combos$avgcvMSE <- apply(cvMSE, 1, mean)
combos$avgcvErrFrac <- apply(cvErrFrac, 1, mean)
## Sqrt units:
## combos$avgsqrtcvMSE <- apply(sqrtcvMSE, 1, mean)
## combos$avgsqrtcvErrFrac <- apply(sqrtcvErrFrac, 1, mean)
## combos$avgorigUnitsqrtcvMSE <- apply(origUnitsqrtcvMSE, 1, mean)
## combos$avgorigUnitsqrtcvErrFrac <- apply(origUnitsqrtcvErrFrac, 1, mean)

write_csv(combos, path="cv_w_degdays_with_rare_at_least_5perc.csv")
## #####################


## #####################
## Make plots of the summary statistics for the various values of
## numBtSampls and numVarSplit.

## Original units:
ggplot(data=combos, aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line() + ggtitle("Cross-validation leaving out groups of degdays")
ggplot(data=combos, aes(x=numBtSamps, y=avgcvErrFrac, color=as.factor(numVarSplit))) + geom_line() + ggtitle("Cross-validation leaving out groups of degdays")
## ggsave("cv_w_degdays_with_rare_at_least_5perc.pdf")

## Sqrt units:
## ggplot(data=combos, aes(x=numBtSamps, y=avgsqrtcvMSE, color=as.factor(numVarSplit))) + geom_line()
## ## X11()
## ggplot(data=combos, aes(x=numBtSamps, y=avgorigUnitsqrtcvMSE, color=as.factor(numVarSplit))) + geom_line()
## ggplot(data=combos, aes(x=numBtSamps, y=avgsqrtcvErrFrac, color=as.factor(numVarSplit))) + geom_line()
## ggplot(data=combos, aes(x=numBtSamps, y=avgorigUnitsqrtcvErrFrac, color=as.factor(numVarSplit))) + geom_line()
## ####################
## ##################################################
