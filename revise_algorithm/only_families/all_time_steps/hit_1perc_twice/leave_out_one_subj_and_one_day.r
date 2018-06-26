library("tidyverse")
library("randomForest")
library("figdim")
library("parallel")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0(taxalevel, "_hit_cutoff_twice_all_time_steps.csv"))
## ##################################################



## ##################################################
## Put the data in wide format; remove days, subj, and rare taxa.

## Move back to wide format.  Unlike our other runs, here we leave
## the subject identifier in the dataset, so that we can try predictions
## with/without designated subjects.
wideT <- taxaT %>%
  filter(taxa!="Rare") %>%
  select(degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay)

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

## Repeated cross-validation runs (1000 of them), leaving out 20% of
## the observations at a time, indicated that the number of variables
## to consider at each split is about 8 (also 9 is very close)
## for the response variable in the original units.
numVarSplit <- 8
## ##################################################



## ##################################################
## Set up these training and valiation datasets, plus matrices to hold
## cross-valiation results.

## We exclude each combination of individual and degree day. This is
## 16 degree days x 6 individuals = 96 combos.
excludeMat <- expand.grid(unique(wideT$subj), unique(wideT$degdays),
                          stringsAsFactors=FALSE)
colnames(excludeMat) <- c("subj", "degdays")
numCVs <- nrow(excludeMat)

## Set up the training and validation datasets corresponding to each
## combo.
crossvalidL <- vector("list", numCVs)
for (i in 1:numCVs){
  lvOut <- (wideT$subj==excludeMat[i,"subj"]) | (wideT$degdays==excludeMat[i,"degdays"])
  trainT <- wideT[!lvOut,] %>% select(-subj)
  validT <- wideT[lvOut,] %>% select(-subj)
  crossvalidL[[i]] <- list(trainT=trainT, validT=validT)
}
rm(i, lvOut, trainT, validT)


## Set up vectors to hold cross-validation results.
cvMSE <- rep(NA, numCVs)
cvErrFrac <- rep(NA, numCVs)
## #########################################



## #########################################
## Set up function for fitting random forest model using original
## units.
origUnitsF <- function(x, mtry, ntree){
  rf <- randomForest(degdays ~ ., data=x$trainT, mtry=mtry, ntree=ntree, importance=T)
  return(predict(rf, newdata=x$validT))
}

## Try using lapply to fit the random forests.
origFitL <- mclapply(crossvalidL, mc.cores=4, origUnitsF, mtry=numVarSplit, ntree=numBtSamps)


## Set up function for fitting random forest model using square root
## units.
## sqrtUnitsF <- function(x, jCombo){
##   sqrtrf <- randomForest(sqrt(degdays) ~ . -subj, data=x$trainT, mtry=combos[jCombo, "numVarSplit"], ntree=combos[jCombo, "numBtSamps"], importance=T)
##   return(predict(sqrtrf, newdata=x$validT))
## }
## #########################################



## #########################################
## Now, calculate the various summary statistics for each cross-validation.

for (i in 1:numCVs){

  ## Get the validation set for this run from the list.
  validT <- crossvalidL[[i]][["validT"]]

  ## Calculate SSTotal for the cross-validation set.
  SSTot <- sum( (validT$degdays-mean(validT$degdays))^2 )
  
  ## Calculate the MSE and error fraction of the SS Total for the
  ## validation data in the original units.
  resid <- origFitL[[i]] - validT$degdays
  cvMSE[i] <- mean(resid^2)
  cvErrFrac[i] <- sum(resid^2)/SSTot
  rm(resid)
}
rm(i, validT, SSTot)
## #########################################



## #########################################
## Make a plot showing how much error we have when we make predictions
## without using a particular day.

residDF <- NULL
for (i in 1:numCVs){

  ## Get the validation set for this run from the list.
  validT <- crossvalidL[[i]][["validT"]]

  ## Calculate the residuals for this validation set.
  resid <- validT$degdays - origFitL[[i]]

  ## Find residuals which correspond to the degree day which was left
  ## out in the training dataset.
  dayOfInterest <- excludeMat[i,"degdays"]
  whichItems <- validT$degdays==dayOfInterest
  residOfInterest <- resid[whichItems]

  ## Build a data frame with these residuals, along with the day and
  ## individual that were left out in this validation.
  iDayDF <- data.frame(degdays=excludeMat[i,"degdays"],
                       subj=excludeMat[i,"subj"],
                       yhat=origFitL[[i]][whichItems],
                       resid=residOfInterest)
  ## Add this data.frame to the what we've already collected.
  residDF <- rbind(residDF, iDayDF)
}
rm(i, validT, resid, dayOfInterest, residOfInterest, iDayDF)


## Write this info out.
write.csv(residDF, file="resids_leave_out_one_subj_and_one_day.csv")
## #########################################




ggplot(residDF, aes(x=degdays, y=resid)) +
  geom_point(aes(col=subj)) +
  labs(x="Degree day", y="Residual")




ggplot(residDF, aes(x=degdays, y=resid)) +
  facet_wrap(~subj) +
  geom_point(aes(col=subj)) +
  labs(x="Actual degree day", y="Residual")

ggplot(residDF, aes(x=yhat, y=resid)) +
  facet_wrap(~subj) +
  geom_point(aes(col=subj)) +
  labs(x="Predicted degree day", y="Residual")
