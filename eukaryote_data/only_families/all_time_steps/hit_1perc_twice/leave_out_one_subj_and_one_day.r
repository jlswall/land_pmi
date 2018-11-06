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

## Repeated cross-validation runs (1000 of them), leaving out 20% of
## the observations at a time, indicated that the number of variables
## to consider at each split is about 15 (also 14 is very close)
## for the response variable in the original units.
numVarSplit <- 15
## ##################################################



## ##################################################
## Set up these training and valiation datasets, plus matrices to hold
## cross-valiation results.

## We exclude each combination of individual and degree day. This is
## 16 degree days x 6 individuals = 96 combos.
excludeCombos <- expand.grid(unique(wideT$subj), unique(wideT$degdays),
                          stringsAsFactors=FALSE)
colnames(excludeCombos) <- c("subj", "degdays")

## With the bacterial data, we had some missing data for certain
## subject-day combos, but we don't have that problem with the
## eukaryotic data.  We do want to turn the matrix above into a tibble.
excludeCombos <- as.tibble(excludeCombos)

## How many times do we want to run each exclusion combo?
numRunsEachCombo <- 1
excludeT <- NULL
for (i in 1:numRunsEachCombo)
    excludeT <- bind_rows(excludeT, excludeCombos)
rm(excludeCombos)
## Order by degdays and then subj.
excludeT <- excludeT %>% arrange(degdays, subj)

## Set up the training and validation datasets corresponding to each
## combo.
numCVs <- nrow(excludeT)
crossvalidL <- vector("list", numCVs)
for (i in 1:numCVs){
  lvOut <- (wideT$subj==pull(excludeT[i,"subj"])) | (wideT$degdays==pull(excludeT[i,"degdays"]))
  trainT <- wideT[!lvOut,] %>% select(-subj)
  validT <- wideT[lvOut,]
  crossvalidL[[i]] <- list(trainT=trainT, validT=validT)
}
rm(i, lvOut, trainT, validT)
## #########################################



## #########################################
## Set up function for fitting random forest model using original
## units.
origUnitsF <- function(x, mtry, ntree){
  rf <- randomForest(degdays ~ ., data=x$trainT, mtry=mtry, ntree=ntree, importance=T)
  return(predict(rf, newdata=x$validT))
}

## Set random seed for reproducibility.
set.seed(416965)

## Try using lapply to fit the random forests.
origFitL <- mclapply(crossvalidL, mc.cores=6, origUnitsF, mtry=numVarSplit, ntree=numBtSamps)
## #########################################



## #########################################
## Collect the residuals, making a note about which day and
## individual were left out.

## Set up vectors to hold cross-validation results.
cvMSE <- rep(NA, numCVs)
cvErrFrac <- rep(NA, numCVs)

residDF <- NULL
for (i in 1:numCVs){

  ## Get the validation set for this run from the list.
  validT <- crossvalidL[[i]][["validT"]]

  ## Calculate SSTotal for the cross-validation set.
  SSTot <- sum( (validT$degdays-mean(validT$degdays))^2 )

  ## Calculate the residuals for this validation set.
  resid <- validT$degdays - origFitL[[i]]

  ## Overall cross-validations statistics, using all residuals.
  cvMSE[i] <- mean(resid^2)
  cvErrFrac[i] <- sum(resid^2)/SSTot
  
  ## ## Find residuals which correspond to the degree day which was left
  ## ## out in the training dataset.
  ## dayOfInterest <- excludeMat[i,"degdays"]
  ## whichItems <- validT$degdays==dayOfInterest
  ## residOfInterest <- resid[whichItems]

  ## Build a data frame with these residuals, along with the day and
  ## individual that were left out in this validation.
  iresidDF <- data.frame(dayOmit=pull(excludeT[i,"degdays"]),
                         subjOmit=pull(excludeT[i,"subj"]),
                         subjactual=validT$subj,
                         yactual=validT$degdays,
                         yhat=origFitL[[i]],
                         resid=resid)
  ## Add this data.frame to the what we've already collected.
  residDF <- rbind(residDF, iresidDF)
}
rm(i, validT, resid, iresidDF)

## Write this info out.
write.csv(residDF, file="resids_leave_out_one_subj_and_one_day.csv", row.names=FALSE)
## #########################################



## #########################################
## We're interested in the error we're likely to get with the model in
## "regular use".  That means that we are interested the prediction
## the model makes for a subject and time slot that we're never
## observed.  That's only ONE prediction of interest per model run.
## So, we look at the root of the mean squared errors for the n
## predictions of interest in these n runs.

myresids <- residDF %>%
  filter((subjactual==subjOmit) & (dayOmit==yactual)) %>%
  pull(resid)
sqrt(mean(myresids^2))
## This is about 243-244, whether I use 1, 10, or 100 runs per
## combo.
## #########################################



## #########################################
## Make plot showing the residuals associated with days which were
## completely left out of the model.

ggplot(residDF %>%
       filter(dayOmit==yactual) %>%
       filter(subjOmit==subjactual),
       aes(x=yactual, y=resid)) +
  geom_point() +
  ## geom_point(aes(col=subjOmit)) +
  geom_hline(yintercept=0) +
  labs(x="Actual degree days", y="Error (actual - estimated)")
ggsave(filename="leave_out_one_subj_and_one_day_residuals.pdf", height=3.5, width=3.5, units="in")

ggplot(residDF, aes(x=yactual, y=resid)) +
  facet_wrap(~subjOmit) +
  geom_point(aes(col=subjOmit)) +
  labs(x="Actual degree day", y="Residual")
## #########################################

