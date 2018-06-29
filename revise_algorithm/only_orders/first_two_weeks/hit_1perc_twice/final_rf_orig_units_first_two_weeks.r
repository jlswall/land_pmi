library("tidyverse")
library("randomForest")
library("figdim")
library("parallel")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "orders"

## Read in cleaned-up phyla, orders, or families taxa.
earlyT <- read_csv(paste0(taxalevel, "_hit_cutoff_twice_first_two_weeks.csv"))
## ##################################################



## ##################################################
## Put the data in wide format and restrict to the first 15 days.

## Move to wide format.
wideT <- earlyT %>%
  filter(taxa!="Rare") %>%
  select(degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay) %>%
  select(-subj)

## Just for reference later, keep the days and degree days, so we can
## look at the time correspondence.
timeT <- earlyT %>% distinct(days, degdays)

rm(earlyT)
## ##################################################



## ##################################################
## From earlier experiments, we figured out these parameters work best
## for the random forest model.

## Number of bootstrap samples.
numBtSamps <- 3000

## Repeated simulations indicated that the number of variables to
## consider at each split is 10 (followed by 11) for response variable
## in the original units.
numVarSplit <- 10
## ##################################################


## ##################################################
## Run the cross-validation for this model, so that we can see what
## the CV MSE looks like.

set.seed(428109)

## Number of times to do cross-validation.
numCVs <- 1000
## How many observations to reserve for testing each time.
numLeaveOut <- round(0.20 * nrow(wideT))


## ###########################
## Set up function for fitting random forest model using original
## units.
origUnitsF <- function(x, mtry, ntree){
  rf <- randomForest(degdays ~ . , data=x$trainT, mtry=mtry, ntree=ntree, importance=T)
  return(predict(rf, newdata=x$validT))
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
origFitL <- mclapply(crossvalidL, mc.cores=4, origUnitsF, mtry=numVarSplit, ntree=numBtSamps)
## ###########################


## ###########################
## For matrix to hold cross-validation results.
cvMSE <- rep(NA, numCVs)
cvErrFrac <- rep(NA, numCVs)

set.seed(743914)


## Now, calculate the various summary statistics for each cross-validation.
residsDF <- NULL
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
  residsDF <- rbind(residsDF, iCaseDF)

  
  ## Calculate the MSE and error fraction of the SS Total for the
  ## validation data in the original units.
  cvMSE[i] <- mean(resid^2)
  cvErrFrac[i] <- sum(resid^2)/SSTot
  rm(resid, iCaseDF)
}
rm(i, validT, SSTot)

write_csv(residsDF, path="final_rf_orig_units_residuals_first_two_weeks.csv")
write_csv(data.frame(cvMSE, cvErrFrac), path="final_rf_orig_units_cvstats_hit_cutoff_twice_first_two_weeks.csv")
rm(cvMSE, cvErrFrac)
## ##################################################



## ##################################################
## Fit the final random forest with all the data (no cross-validation).

set.seed(5580532)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(degdays ~ . , data=wideT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

init.fig.dimen(file=paste0("orig_units_hit_cutoff_twice_first_two_weeks_orders_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of order taxa (orig. units, first 2 weeks)")
dev.off()


## Find residuals:
resids <- rf$predicted - wideT$degdays

## Print out RMSE:
sqrt( mean( resids^2 ) )
## RMSE: 54.67932

## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (wideT$degdays - mean(wideT$degdays))^2 ) )
## Expl. frac.: 0.8724671
## ##################################################



## ##################################################
## Make graph of just IncNodePurity alone.

## Turn importance measures into a tibble, sorted by IncNodePurity in
## increasing order.
importanceT <- importance(rf) %>%
  as.data.frame() %>% as_tibble() %>%
  rownames_to_column("order") %>%
  arrange(IncNodePurity)
## Turn order names into factors, so that we can make the bar chart
## with the bars in decreasing order.
importanceT$order <- factor(importanceT$order, levels=importanceT$order)
ggplot(importanceT %>% top_n(10, wt=IncNodePurity),
       aes(x=order, y=IncNodePurity)) +
  coord_flip() +
  geom_col() +
  labs(x="Order", y="Decrease in node impurity")
ggsave(filename="orig_units_first_two_weeks_orders_barchart.pdf", height=2.5, width=4, units="in")
## ##################################################
