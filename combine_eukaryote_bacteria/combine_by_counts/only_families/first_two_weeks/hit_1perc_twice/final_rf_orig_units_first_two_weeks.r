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

## Move back to wide format.
wideT <- taxaT %>%
  filter(taxa!="Rare") %>%
  select(degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay) %>%
  select(-subj)

## Just for reference later, keep the days and degree days, so we can
## look at the time correspondence.
timeT <- taxaT %>% distinct(days, degdays)

## rm(taxaT)  ## Use to make plot of influential taxa at finish.
## ##################################################



## ##################################################
## From earlier experiments, we figured out these parameters work best
## for the random forest model.

## Number of bootstrap samples.
numBtSamps <- 3000

## Repeated cross-validation runs, leaving out 20% of the observations
## at a time, indicated that the number of variables to consider at
## each split is about 15 (with 14, 16, and 13 close) for the response
## variable in the original units.
numVarSplit <- 15
## ##################################################



## ##################################################
## Run the cross-validation for this model, so that we can see what
## the CV MSE looks like.

set.seed(6091942)

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

set.seed(7147038)


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

  
  ## Calculate the MSE and error fraction of the SS Total for the
  ## validation data in the original units.
  cvMSE[i] <- mean(resid^2)
  cvErrFrac[i] <- sum(resid^2)/SSTot
  rm(resid, iCaseDF)
}
rm(i, validT, SSTot)

write_csv(residDF, path="final_rf_orig_units_residuals_first_two_weeks.csv")
write_csv(data.frame(cvMSE, cvErrFrac), path="final_rf_orig_units_cvstats_first_two_weeks.csv")
rm(cvMSE, cvErrFrac)
## ##################################################



## ##################################################
## Fit the final random forest with all the data (no cross-validation).

set.seed(258555)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(degdays ~ . , data=wideT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

init.fig.dimen(file=paste0("orig_units_first_two_weeks_families_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of family taxa (orig. units, first two weeks)")
dev.off()


## Find residuals:
resids <- rf$predicted - wideT$degdays

## Print out RMSE:
sqrt( mean( resids^2 ) )
## RMSE: 57.19763

## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (wideT$degdays - mean(wideT$degdays))^2 ) )
## Expl. frac.: 0.8604493
## ##################################################



## ##################################################
## Make graph of just IncNodePurity alone.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 8

## Turn importance measures into a tibble, sorted by IncNodePurity in
## increasing order.
importanceT <- importance(rf) %>%
  as.data.frame() %>% as_tibble() %>%
  rownames_to_column("family") %>%
  arrange(IncNodePurity)
## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
importanceT$family <- factor(importanceT$family, levels=importanceT$family)
ggplot(importanceT %>% top_n(n, wt=IncNodePurity),
       aes(x=family, y=IncNodePurity)) +
  coord_flip() +
  geom_col() +
  labs(x="Bac. and euk. family-level taxa", y="Decrease in node impurity") +
  theme(axis.title=element_text(size=rel(0.8)), axis.text=element_text(size=rel(0.8))) 
ggsave(filename="orig_units_first_two_weeks_families_barchart.pdf", height=2.5, width=4.5, units="in")
## ##################################################



## ##################################################
## Make graph of just %IncMSE alone.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 8

## Turn importance measures into a tibble, sorted by IncNodePurity in
## increasing order.
importanceT <- importance(rf) %>%
  as.data.frame() %>% as_tibble() %>%
  rownames_to_column("family") %>%
  arrange(`%IncMSE`)
## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
importanceT$family <- factor(importanceT$family, levels=importanceT$family)
ggplot(importanceT %>% top_n(n, wt=`%IncMSE`),
       aes(x=family, y=`%IncMSE`)) +
  coord_flip() +
  geom_col() +
  labs(x="Bac. and euk. family-level taxa", y="Mean % decrease in MSE when excluded") +
  theme(axis.title=element_text(size=rel(0.8)), axis.text=element_text(size=rel(0.8))) 
ggsave(filename="orig_units_first_two_weeks_families_PercIncMSE_barchart.pdf", height=2.5, width=4.5, units="in")
## ##################################################



## ##################################################
## Make scatter plots of the percentages over time (by taxa) for the
## top n taxa in terms of %IncMSE.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 8

## Save the names of the families that are in the top 10 in
## terms of %IncMSE.
topChoices <- as.character(importanceT %>% arrange(desc(`%IncMSE`)) %>% pull(family))[1:n]

## Find the percentages for these taxa.
chooseT <- taxaT %>%
  filter(taxa %in% topChoices)
chooseT$taxa <- factor(chooseT$taxa, levels=topChoices)

ggplot(chooseT, aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  labs(x="Degree days", y="Fraction", color="Cadaver") +
  theme(legend.title=element_text(size=rel(0.8)), legend.text=element_text(size=rel(0.8))) + 
  ## Allow diff. y-scales across panels.
  facet_wrap(~taxa, ncol=4, scales="free_y") 
  ## facet_wrap(~taxa)  ## Keep y-scales same across panels.
ggsave("infl_comb_family_first_two_weeks_scatter.pdf", width=8, height=4, units="in")
## ##################################################



## ##################################################
## Make plot of residuals.

ggplot(residDF, aes(x=yactual, y=resid)) +
  geom_point() +
  geom_hline(yintercept=0) + 
  labs(x="Actual accumulated degree days", y="Error (actual - estimated)")
ggsave(filename="orig_units_first_two_weeks_families_residuals.pdf", height=3.5, width=3.5, units="in")
## ##################################################
