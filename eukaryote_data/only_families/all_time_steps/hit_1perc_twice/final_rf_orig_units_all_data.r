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

## Move back to wide format.
wideT <- taxaT %>%
  filter(taxon!="Rare") %>%
  select(degdays, subj, taxon, fracBySubjDay) %>%
  spread(taxon, fracBySubjDay) %>%
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

## Repeated cross-validation runs (1000 of them), leaving out 20% of
## the observations at a time, indicated that the number of variables
## to consider at each split is about 15 (also 14 is very close)
## for the response variable in the original units.
numVarSplit <- 15
## ##################################################



## ##################################################
## Fit the random forest using all the data (do not hold any back for
## cross-validation).

## ###########################
## Set up function for fitting random forest model using full dataset.
fitFullDataF <- function(x, mtry, ntree){

  rf <- randomForest(degdays ~ . , data=x, mtry=mtry, ntree=ntree, importance=T)

  ## Order the taxa according to decreasing values of the scaled
  ## importance %IncMSE.  Return as a tibble, with %IncMSE column
  ## renamed to PercIncMSE.
  importanceT <- as.tibble(importance(rf), rownames="taxa") %>% rename(PercIncMSE=`%IncMSE`) %>% arrange(desc(PercIncMSE))

  myresults <- list(predicted=rf[["predicted"]], y=rf[["y"]], importanceT=importanceT)

  return(myresults)
}
## ###########################


## ###########################
## Set up list as input for the function, which is just the original
## dataset repeated over and over.
numRepeat <- 1000
repDataL <- vector("list", numRepeat)
for (i in 1:numRepeat)
  repDataL[[i]] <- wideT
## ###########################


## ###########################
## Now, fit random forests to the full dataset over and over.

set.seed(718460)
fullResultsL <- mclapply(repDataL, mc.cores=6, fitFullDataF, mtry=numVarSplit, ntree=numBtSamps)

## Calculate the RMSE and pseudo-Rsquared for these runs with the full
## dataset.
fullRMSE <- rep(NA, 1000)
fullRsq <- rep(NA, 1000)
fullImportanceT <- NULL
for (i in 1:length(fullResultsL)){

  ## Get results from run i.
  iTmp <- fullResultsL[[i]]

  ## Find residuals:
  resids <- iTmp$predicted - iTmp$y

  ## Calculate the RMSE and pseudo-Rsquared:
  fullRMSE[i] <- sqrt( mean( resids^2 ) )
  fullRsq[i] <- 1.0 - ( sum(resids^2)/sum( (iTmp$y - mean(iTmp$y))^2 ) )

  ## Store measures of importance in a long tibble.
  fullImportanceT <- rbind(fullImportanceT, iTmp$importanceT)
}
rm(i, iTmp, resids)

## See summary of %IncMSE (measure of importance) over all model runs.
fullImportanceT %>% group_by(taxa) %>% summarize(meanPercIncMSE=mean(PercIncMSE), lbPercIncMSE=quantile(PercIncMSE, 0.025), ubPercIncMSE=quantile(PercIncMSE, 0.975)) %>% arrange(desc(meanPercIncMSE))

## Get summary statistics for report.
c(mean(fullRMSE), 1.96*sd(fullRMSE))
## RMSE: 177.680949   2.411631
c(mean(fullRsq), 1.96*sd(fullRsq))
## Rsq: 0.895643282 0.002831568

write_csv(data.frame(fullRMSE, fullRsq), path="cvstats_w_full_dataset_final_params.csv")
rm(fullRMSE, fullRsq)
## ###########################
## ##################################################



## ##################################################
## Fit the final random forest with all the data (no cross-validation).

set.seed(2240635)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(degdays ~ . , data=wideT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)

## init.fig.dimen(file=paste0("orig_units_all_data_families_imp_plot.pdf"), width=8, height=6)
## varImpPlot(rf, main="Importance of eukaryotic family-level taxa (all time steps)")
## dev.off()


## Find residuals:
resids <- rf$predicted - wideT$degdays

## Print out RMSE:
sqrt( mean( resids^2 ) )
## RMSE: 177.5454

## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (wideT$degdays - mean(wideT$degdays))^2 ) )
## Expl. frac.: 0.8958074
## ##################################################



## ## ##################################################
## ## Make graph of just IncNodePurity alone.

## ## Get the top "n" (whether 8, 10, whatever) influential taxa.
## n <- 8

## ## Turn importance measures into a tibble, sorted by IncNodePurity in
## ## increasing order.
## importanceT <- importance(rf) %>%
##   as.data.frame() %>% as_tibble() %>%
##   rownames_to_column("family") %>%
##   arrange(IncNodePurity)
## ## Turn family names into factors, so that we can make the bar chart
## ## with the bars in decreasing order.
## importanceT$family <- factor(importanceT$family, levels=importanceT$family)
## ggplot(importanceT %>% top_n(n, wt=IncNodePurity),
##        aes(x=family, y=IncNodePurity)) +
##   coord_flip() +
##   geom_col() +
##   labs(x="Eukaryotic family-level taxa", y="Decrease in node impurity")
## ggsave(filename="orig_units_all_data_families_IncNodePurity_barchart.pdf", height=2.5, width=4.5, units="in")
## ## ##################################################



## ##################################################
## Make graph of just %IncMSE alone.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 6

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
  labs(x="Eukaryotic family-level taxa", y="Mean % decrease in MSE when excluded")
ggsave(filename="orig_units_all_data_families_PercIncMSE_barchart.pdf", height=2.5, width=4.5, units="in")
## ##################################################



## ##################################################
## Make scatter plots of the percentages over time (by taxa) for the
## top n taxa in terms of %IncMSE.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 6

## Save the names of the families that are in the top 10 in
## terms of %IncMSE.
topChoices <- as.character(importanceT %>% arrange(desc(`%IncMSE`)) %>% pull(family))[1:n]

## Find the percentages for these taxa.
chooseT <- taxaT %>%
  filter(taxon %in% topChoices)
chooseT$taxon <- factor(chooseT$taxon, levels=topChoices)

ggplot(chooseT, aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  labs(x="Degree days", y="Fraction", color="Cadaver") +
  theme(legend.title=element_text(size=rel(0.8)), legend.text=element_text(size=rel(0.8))) + 
  ## Allow diff. y-scales across panels.
  facet_wrap(~taxon, ncol=3, scales="free_y") 
  ## facet_wrap(~taxa)  ## Keep y-scales same across panels.
ggsave("infl_euk_family_all_data_scatter.pdf", width=8, height=4, units="in")
## ##################################################



## ## ##################################################
## ## Make plot of residuals.

## ggplot(residDF, aes(x=yactual, y=resid)) +
##   geom_point() +
##   geom_hline(yintercept=0) + 
##   labs(x="Actual accumulated degree days", y="Error (actual - estimated)")
## ggsave(filename="orig_units_all_data_families_residuals.pdf", height=3.5, width=3.5, units="in")
## ## ##################################################




## ##################################################
## There is a break between the top 5 and the rest in terms of their
## influence, as measured by %IncMSE or IncNodePurity.

## Get these top 5 influential taxa.
n <- 5

## Turn importance measures into a tibble, sorted by IncNodePurity in
## increasing order.
importanceT <- importance(rf) %>%
  as.data.frame() %>% as_tibble() %>%
  rownames_to_column("family") %>%
  arrange(`%IncMSE`)


## Save the names of the families that are in the top 10 in
## terms of %IncMSE.
topChoices <- as.character(importanceT %>% arrange(desc(`%IncMSE`)) %>% pull(family))[1:n]

## Find the percentages for these taxa.
chooseT <- taxaT %>%
  filter(taxon %in% topChoices)

## Average the value across cadavers for each taxa and each day.
summaryChoiceT <- chooseT %>% group_by(taxon, days, degdays) %>% summarize(meanPercByDay=100*mean(fracBySubjDay), medianPercByDay=100*median(fracBySubjDay))

ggplot(summaryChoiceT, aes(x=degdays, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon)) +
  labs(x="Accumulated degree days (ADD)", y="Average percentage composition", title="Family-level, all time steps", tag="A") +
  scale_x_continuous(breaks=unique(summaryChoiceT$degdays)) + 
  theme_bw()

ggplot(summaryChoiceT, aes(x=degdays, medianFracByDay)) +
  geom_line(aes(color=taxon)) +
  theme_bw()
## ##################################################
