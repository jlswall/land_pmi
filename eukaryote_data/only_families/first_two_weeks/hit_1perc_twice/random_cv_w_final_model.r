library("tidyverse")
library("randomForest")
library("figdim")
library("parallel")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0(taxalevel, "_hit_cutoff_twice_first_two_weeks.csv"))
## ##################################################



## ##################################################
## Put the data in wide format; remove days, subj, and rare taxa.

## Move back to wide format.  Unlike our other runs, here we leave the
## subject identifier in the dataset, so that we can try predictions
## with/without designated subjects.
wideT <- taxaT %>%
  filter(taxon!="Rare") %>%
  select(degdays, subj, taxon, fracBySubjDay) %>%
  spread(taxon, fracBySubjDay)

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

## Repeated cross-validation runs, leaving out 20% of
## the observations at a time, indicated that the number of variables
## to consider at each split is about 15 (with 14 and 16 close)
## for the response variable in the original units.
numVarSplit <- 15
## ##################################################



## ##################################################
## Run the cross-validation for this model, so that we can see what
## the CV MSE looks like.

## Number of times to do cross-validation.
numCVs <- 1000

## How many observations to reserve for testing each time.  This
## number is chosen to be same as the number that are used in leave
## out one subject and day cross-validation.  That is, we leave out chosen
## subject on each each of 10 days, leave out the 5 other individuals
## on the chosen day.
numLeaveOut <- 15


## ###########################
## Set up function for fitting random forest model using original
## units for cross-validation.
origUnitsF <- function(x, mtry, ntree){
  rf <- randomForest(degdays ~ . , data=x$trainT, mtry=mtry, ntree=ntree, importance=T)
  return(predict(rf, newdata=x$validT))
}
## ###########################


## ###########################
## Get set up for cross-validation.
set.seed(9820929)
crossvalidL <- vector("list", numCVs)
for (i in 1:numCVs){
  lvOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)
  trainT <- wideT[-lvOut,] %>% select(-subj)
  validT <- wideT[lvOut,] %>% select(-subj)
  crossvalidL[[i]] <- list(trainT=trainT, validT=validT)
}
rm(i, lvOut, trainT, validT)
## ###########################


## ###########################
## Conduct cross-validation.
set.seed(957703)
origFitL <- mclapply(crossvalidL, mc.cores=6, origUnitsF, mtry=numVarSplit, ntree=numBtSamps)
## ###########################


## ###########################
## Find MSE and pseudo-Rsquared for cross-validation results.
cvRMSE <- rep(NA, numCVs)
cvRsq <- rep(NA, numCVs)

## Now, calculate the various summary statistics for each cross-validation.
residDF <- NULL
for (i in 1:numCVs){
  ## Get the validation set for this run from the list.
  validT <- crossvalidL[[i]][["validT"]]

  ## Calculate SSTotal for the cross-validation set.
  SSTot <- sum( (validT$degdays-mean(validT$degdays))^2 )

  ## Calculate the residuals for this validation set.
  resid <- validT$degdays - origFitL[[i]]

  ## Build a data frame with the actual response and the estimated
  ## response.
  iCaseDF <- data.frame(yactual=validT$degdays, yhat=origFitL[[i]],
                        resid=resid)
  ## Add this data frame to what we've already collected.
  residDF <- rbind(residDF, iCaseDF)
  
  ## Calculate the RMSE and pseudo-Rsquared for the validation data
  ## (in the original units).
  cvRMSE[i] <- sqrt( mean(resid^2) )
  cvRsq[i] <- 1.0 - ( sum(resid^2)/SSTot )
  rm(resid, iCaseDF)
}
rm(i, validT, SSTot)

## Get summary statistics for report.
c(mean(cvRMSE), 1.96*sd(cvRMSE))
## RMSE: 63.86031 +/- 26.41300

## write_csv(residDF, path="final_rf_orig_units_residuals_first_two_weeks.csv")
write_csv(data.frame(cvRMSE, cvRsq), path="cvstats_random_w_final_params.csv")
rm(cvRMSE, cvRsq)
## ##################################################



## ##################################################
## Make plot of residuals.

ggplot(residDF, aes(x=yactual, y=resid)) +
  geom_point() +
  geom_hline(yintercept=0) + 
  labs(x="Actual accumulated degree days", y="Error (actual - estimated)")
ggsave(filename="resids_random_cv_w_final_params.pdf", height=3.5, width=3.5, units="in")
## ##################################################
