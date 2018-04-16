library("tidyverse")
library("randomForest")
library("readr")

## Load the random forest model that we fit using the first 2 weeks of
## the original data.
load("../revise_algorithm/only_orders/first_two_weeks/final_rf_orders_first_two_weeks.RData")

## Get list of variables considered in the model.
rfCovar <- attr(rf[["terms"]], "term.labels")


## Read in the validation data.
validT <- read_csv("orders_massaged.csv")


## ##################################################
## Compare the percentages with those from the original dataset for
## the most important variables.

## Find the variables by order of importance in original random forest
## model.
orderIncMSE <- order(rf[["importance"]][,"%IncMSE"], decreasing=TRUE)
varsByImportance <- rownames(rf[["importance"]][orderIncMSE,])
rm(orderIncMSE)

## Read in the original data.
origT <- read_csv("../revise_algorithm/orders_massaged.csv")


## Put on dots to represent validation measurements.
mypch <- c(1, 15, 16)
names(mypch) <- unique(validT$subj)
mycolors <- c("blue", "darkgreen", "slategray")
names(mycolors) <- unique(validT$subj)

pdf(file="compare_taxa_orig_vs_valid_data.pdf")
par(mfrow=c(2,1))

for (iTaxa in varsByImportance){

  ## Subset and draw max/min lines for original data.
  subT <- origT %>%
    filter((taxa==iTaxa) & (degdays<=120)) %>%
    group_by(degdays) %>%
    summarize(minFrac=min(fracBySubjDay), maxFrac=max(fracBySubjDay))

  ## Subset the validation measurements.
  subValidT <- validT %>%
    filter(taxa==iTaxa)

  ## Build graph.
  ## First, find range of counts.
  fracRng <- range(c(subT[,c("minFrac", "maxFrac")], subValidT[,"fracBySubjDay"]))
  plot(range(c(subT[,"degdays"], subValidT[,"degdays"])), fracRng,
       type="n", mar=c(3, 3, 2, 1), xlab="Degree days",
       ylab="Fractional composition")
  title(iTaxa, line=1)
  for (i in 1:nrow(subT))
    with(subT[i,], lines(rep(degdays, 2), c(minFrac, maxFrac)))  
  for (iSubj in unique(subValidT$subj))
    with(subValidT %>% filter((subj==iSubj) & (taxa==iTaxa)),
         points(degdays, fracBySubjDay,
                pch=mypch[iSubj], col=mycolors[iSubj])
         )
  legend("topright", names(mypch), pch=mypch, col=mycolors, cex=0.5, title="Valid. pigs")
}
dev.off()
## ##################################################




## For the validation data, assess taxa variability across the first
## few days, averaged across pigs, in terms of fractions.
avgTaxaOverdegdaysT <- validT %>%
  filter(degdays <= 150) %>%
  group_by(degdays, taxa) %>%
  summarize(avgFracByDegday=mean(fracBySubjDay)) %>%
  filter(avgFracByDegday >= 0.01)



## ##################################################
## Put the data in wide format.

## Move back to wide format.
wideT <- validT %>%
  select(degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay) %>%
  select(-subj)


## There's one taxa from the original model which isn't in the
## validation data.  That is "Streptophyta".
rfCovar[!(rfCovar %in% colnames(wideT))]
## Add a column of this name with zeroes in it.
wideT$Streptophyta <- 0.0
## ##################################################


## Make prediction for the validation set, using the original model.
predictValid <- predict(rf, newdata=wideT)
## Find residuals:
resids <- predictValid - wideT$degdays
## Print out RMSE:
sqrt( mean( resids^2 ) )
## Estimate of explained variance, which R documentation calls "pseudo
## R-squared"
1 - ( sum(resids^2)/sum( (wideT$degdays - mean(wideT$degdays))^2 ) )

## The predictions are not good, and I think this is because the
## fractions for the various taxa are quite different.


