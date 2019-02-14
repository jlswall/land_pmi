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
## Fit the final random forest with all the data (no cross-validation).

set.seed(2240635)

## Fit the random forest model on all the data (no cross-validation).
rf <- randomForest(degdays ~ . , data=wideT, mtry=numVarSplit,
                   ntree=numBtSamps, importance=T)


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




## ##################################################
## Make four-panel figure for use in publication.


## ########################
## Make graph of just %IncMSE alone.

## Get the top "n" (whether 8, 10, whatever) influential taxa.
n <- 10

## Turn importance measures into a tibble, sorted by IncNodePurity in
## increasing order.
importanceT <- importance(rf) %>%
  as.data.frame() %>% as_tibble() %>%
  rownames_to_column("family") %>%
  arrange(`%IncMSE`)
## Turn family names into factors, so that we can make the bar chart
## with the bars in decreasing order.
importanceT$family <- factor(importanceT$family, levels=importanceT$family)
tlPanel <- ggplot(importanceT %>% top_n(n, wt=`%IncMSE`),
                  aes(x=family, y=`%IncMSE`)) +
  theme_minimal() +
  coord_flip() +
  geom_col() +
  labs(x="Eukaryotic Family-level Taxa", y="Mean % Decrease in MSE When Taxa Excluded")
## ########################


## ########################
## Show line plot of changing relative abundance for the top 5 taxa.

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
summTopT <- chooseT %>% group_by(taxon, days, degdays) %>% summarize(meanPercByDay=100*mean(fracBySubjDay), medianPercByDay=100*median(fracBySubjDay))

## In Luisa's paper, she had a plot of average relative abundance
## vs. time for 5 taxa that she identified as being present in large
## numbers.  Her x-axis had the time steps evenly spaced (not
## reflecting actual time passage), with each tick mark labeled with
## the day/degreeday.  To do this, but yet keep days in order, we need
## to build a new factor variable of the form day/degree day, with
## ordered levels.
orderedLevels <- with(timeT, paste(degdays, days, sep="/"))
summTopT$dayADD <- factor(with(summTopT, paste(degdays, days, sep="/")), levels=orderedLevels)
rm(orderedLevels)
dev.new(width=4.5, height=4)
trPanel <- ggplot(summTopT, aes(x=dayADD, y=meanPercByDay, group=taxon)) +
  geom_line(size=1.25, aes(color=taxon)) +
  scale_y_continuous(limits=c(0, 100), expand=c(0,0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5),
        legend.position=c(0.95, 0.98),
        legend.justification=c("right", "top"),
        legend.title=element_blank(),
        legend.key.size=unit(0.5, 'lines'),
        legend.background=element_rect(fill="white")) +
  labs(x="Accumulated Degree Days/Days", y="Relative Abundance")## tag="A")
## ########################


## ########################
## This panel is based on Fig. 1c from Pechal et al. (2015).  It
## is predicted vs. actual ADD.

## Make a tibble of actual and predicted values for each observation.
predvactT <- as.tibble(data.frame(predicted=rf$predicted, actual=rf$y))
Rsq <- with(predvactT, round(cor(actual, predicted)^2, 2))
## RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(resids^2)), 2)  
blPanel <- ggplot(predvactT, aes(x=actual, y=predicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=50, y=1700, hjust=0, label=paste("R^2  ==", Rsq), parse=T) +
  annotate("text", x=50, y=1600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(0, max(as.vector(predvactT))), y=c(0, max(as.vector(predvactT)))) +
  labs(x="Actual Accumulated Degree Days", y="Predicted Accumulated Degree Days")
## ########################


## ########################
## This panel shows the same results as the panel above, but on the
## log scale.

## Make new columns with log, base 10.  We have 6 observations made at
## day/ADD 0/0.  Since log10(0) is undefined, we make these values NA.
predvactT$logactual <- with(predvactT, ifelse(actual>0, log10(actual), NA))
predvactT$logpredicted <- with(predvactT, ifelse(predicted>0, log10(predicted), NA))
minAxisLmt <- min(c(predvactT$logactual, predvactT$logpredicted), na.rm=T)
maxAxisLmt <- max(c(predvactT$logactual, predvactT$logpredicted), na.rm=T)

Rsq <- with(predvactT, round(cor(actual, predicted)^2, 2))
## RMSE around 1:1 line, not regression line.
RMSE <- round(sqrt(mean(resids^2)), 2)  
brPanel <- ggplot(predvactT, aes(x=logactual, y=logpredicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  ## annotate("text", x=50, y=1700, hjust=0, label=paste("R^2  ==", Rsq), parse=T) +
  ## annotate("text", x=50, y=1600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(minAxisLmt, maxAxisLmt), y=c(minAxisLmt, maxAxisLmt)) +
  labs(x="Log 10 of Actual Accumulated Degree Days", y="Log 10 of Predicted Accumulated Degree Days")
## ########################


## ########################
library("cowplot")
## Try theme(plot.margin) when creating graphs to add some space.
plot_grid(tlPanel, trPanel, blPanel, brPanel, labels=c("a", "b", "c", "d"))##, rel_widths=c(1.125, 1), rel_heights=c(1, 1))
ggsave(file="relative_abundance_Rsq_rmse.pdf", width=8.5, height=8.5, units="in")
## ########################
## ##################################################



## ##################################################
## Tal also wanted to see the predicted vs. actual scatterplot in log mode.

## Make new columns with natural log.  For values that are 0, the log
## is undefined.  I make these values 0.
predvactT$logactual <- with(predvactT, ifelse(actual>0, log(actual), 0))
predvactT$logpredicted <- with(predvactT, ifelse(predicted>0, log(predicted), 0))
minAxisLmt <- min(c(predvactT$logactual, predvactT$logpredicted), na.rm=T)
maxAxisLmt <- max(c(predvactT$logactual, predvactT$logpredicted), na.rm=T)
Rsq <- with(predvactT, round(cor(logactual, logpredicted)^2, 2))
ggplot(predvactT, aes(x=logactual, y=logpredicted)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=0.5, y=6.5, hjust=0, label=paste("R^2  ==", Rsq), parse=T) +
##  annotate("text", x=50, y=1600, hjust=0, label=paste("RMSE = ", RMSE)) + 
  coord_fixed(ratio=1) +
  theme_bw() + 
  lims(x=c(minAxisLmt, maxAxisLmt), y=c(minAxisLmt, maxAxisLmt)) +
  labs(x="Natural log of actual accumulated degree days", y="Natural log of predicted accumulated degree days")##tag="B")
ggsave(file="scatterplot_log_actual_predicted.pdf", width=4, height=4, units="in")
