library("tidyverse")
library("randomForest")
library("figdim")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0(taxalevel, "_massaged.csv"))
## ##################################################



## ##################################################
## Put the data in wide format.

## Move back to wide format.
wideT <- taxaT %>%
  ungroup() %>%
  select(days, degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay)
## ##################################################



## ##################################################
## Try using all the data, even though the observations at the end are
## taken less often than the ones at the beginning.

## ##########################
## Look at a cluster analysis, just to see whether observations made
## at the same time tend to be grouped together.  Results are mixed.
set.seed(4517802)
allT <- wideT %>% select(-subj, -Rare, -days)
hc.out <- hclust(dist(allT %>% select(-degdays)), method="average")
plot(hc.out, labels=allT$degdays)
## ##########################


## ##########################
## Try random forests for regression using "days" as the response
## variable.

## #########
## How many predictors?  (All columns except response: "degdays").
numPredictors <- ncol(allT) - 1

## Try different numbers of bootstrap samples.
numBtSampsVec <- seq(200, 1000, by=200)
## numBtSampsVec <- seq(500, 4000, by=500)

## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
## numVarSplitVec <- sort(unique(c(floor(sqrt(numPredictors)), ceiling(sqrt(numPredictors)), floor(numPredictors/3), ceiling(numPredictors/3), floor(numPredictors/2), numPredictors)))
numVarSplitVec <- c(floor(sqrt(numPredictors)):ceiling(numPredictors/3))

## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)

## ############
## ## Uncomment the following to do cross-validation over and over.
## ## Number of times to do cross-validation.
## numCVs <- 100
## ## How many observations to reserve for testing each time.
## numLeaveOut <- round(0.10 * nrow(wideT))
## ## For matrix to hold cross-validation results.
## cvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
## ############

## Number of times to do cross-validation.
numFolds <- 10
## Figure out which observations are in which fold by randomly
## assigning the numbers 1-10 to the various rows.  There are 93
## observations, so we assign an extra 1, 2, 3.
set.seed(527689)
numbersToAssign <- c( rep( 1:10, floor(nrow(allT)/10) ), 1:(nrow(allT) %% 10) )
whichFold <- sample(numbersToAssign, replace=F)

## For matrix to hold cross-validation results.
cvMSE <- matrix(NA, nrow(combos), ncol=numFolds)
cvSSM <- matrix(NA, nrow(combos), ncol=numFolds)
sqrtcvMSE <- matrix(NA, nrow(combos), ncol=numFolds)
sqrtcvSSM <- matrix(NA, nrow(combos), ncol=numFolds)


## Do cross-validation.
for (i in 1:numFolds){
  ## ############
  ## for (i in 1:numCVs){
  ## ##  Uncomment the following to do cross-validation over and over.
  ## ## Determine training and cross-validation set.
  ## whichLeaveOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)    
  ## subT <- wideT[-whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
  ## cvsetT <- wideT[whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
  ## ############

  ## Leave out all the rows assigned to i.
  whichLeaveOut <- whichFold==i
  cvsetT <- allT[whichLeaveOut,]
  subT <- allT[!whichLeaveOut,]
  
  for (j in 1:nrow(combos)){    
    rf <- randomForest(degdays ~ . , data=subT, mtry=combos[j,"numVarSplit"], ntree=combos[j,"numBtSamps"], importance=T)
    fitTest <- predict(rf, newdata=cvsetT)
    cvMSE[j,i] <- mean((fitTest - cvsetT$degdays)^2)
    cvSSM[j,i] <- sum((fitTest - cvsetT$degdays)^2)/sum((cvsetT$degdays - mean(cvsetT$degdays))^2)
  }

  ## In sqrt units:
  for (j in 1:nrow(combos)){    
    sqrtrf <- randomForest(sqrt(degdays) ~ . , data=subT, mtry=combos[j,"numVarSplit"], ntree=combos[j,"numBtSamps"], importance=T)
    fitTest <- predict(sqrtrf, newdata=cvsetT)
    sqrtcvMSE[j,i] <- mean((fitTest^2 - cvsetT$degdays)^2)
    sqrtcvSSM[j,i] <- sum((fitTest^2 - cvsetT$degdays)^2)/sum((cvsetT$degdays - mean(cvsetT$degdays))^2)
  }
}
rm(i, j, fitTest, rf, sqrtrf)

combos$avgcvMSE <- apply(cvMSE, 1, mean)
combos$avgsqrtcvMSE <- apply(sqrtcvMSE, 1, mean)
combos$avgcvSSM <- apply(cvSSM, 1, mean)
combos$avgsqrtcvSSM <- apply(sqrtcvSSM, 1, mean)
write.csv(combos, file="all_data_avg_cv_metrics.csv")
## combos$q1cvMSE <- apply(cvMSE, 1, quantile, 0.25)
## combos$q3cvMSE <- apply(cvMSE, 1, quantile, 0.75)
ggplot(data=combos, aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) +
  geom_line()
X11()
ggplot(data=combos, aes(x=numBtSamps, y=avgsqrtcvMSE, color=as.factor(numVarSplit))) + geom_line()
## plot(range(numBtSampsVec),
##      range(as.vector(combos[,c("q1cvMSE", "q3cvMSE")])), type="n")
## my.colors <- c("black", "red", "green", "blue", "orange", "cyan", "magenta", "yellow", "brown")
## for (i in c(1:length(numVarSplitVec))){
##   subDF <- subset(combos, numVarSplit==numVarSplitVec[i])
##   with(subDF, lines(numBtSamps, avgcvMSE, col=my.colors[i]))
##   with(subDF, lines(numBtSamps, q1cvMSE, lty=2, col=my.colors[i]))
##   with(subDF, lines(numBtSamps, q3cvMSE, lty=2, col=my.colors[i]))  
## }
## rm(subDF, i)

yrng <- range(combos[,"avgcvMSE"])
ggplot(data=subset(combos, numVarSplit<=5), aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line() + ylim(yrng)
X11()
ggplot(data=subset(combos, (numVarSplit>5) & (numVarSplit<=8)), aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line() + ylim(yrng)
X11()
ggplot(data=subset(combos, numVarSplit>8), aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line() + ylim(yrng)
## For orders: mtry=6-8 and ntrees >= 3000
## For families: mtry=4 and ntrees >= 3000
## For phyla: mtry=2 and ntrees >= 2500
## ####################


## ####################
## Given the above, use a setting of mtry=8 and ntrees=5000 to get an
## typical answer we can use in graphics.

## #####
## THIS IS FOR ORDER-LEVEL TAXA.
set.seed(5434890)
rf <- randomForest(days ~ . -subj -Rare -degdays, data=wideT, mtry=8, ntree=3000, importance=T)
init.fig.dimen(file=paste0(taxalevel, "_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of order-level taxa")
dev.off()
## #####


## #####
## THIS IS FOR FAMILY-LEVEL TAXA.
set.seed(5434890)
rf <- randomForest(days ~ . -subj -Rare -degdays, data=wideT, mtry=4, ntree=3000, importance=T)
init.fig.dimen(file=paste0(taxalevel, "_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of family-level taxa")
dev.off()
## #####


## #####
## THIS IS FOR PHYLUM-LEVEL TAXA.
set.seed(4804790)
rf <- randomForest(days ~ . -subj -Rare -degdays, data=wideT, mtry=2, ntree=3000, importance=T)
init.fig.dimen(file=paste0(taxalevel, "_imp_plot.pdf"), width=8, height=6)
varImpPlot(rf, main="Importance of phylum-level taxa")
dev.off()
## #####
## ##################################################


## ##################################################
## Try random forests for regression using "degdays" as the response
## variable.

## How many predictors are there?  (All taxa except "rare".)
numPredictors <- nrow(taxaT %>% filter(taxa!="Rare") %>% distinct(taxa))

## Try different numbers of bootstrap samples.
numBtSampsVec <- seq(500, 4000, by=500)
## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
## numVarSplitVec <- 6:8
## numVarSplitVec <- unique(c(floor(sqrt(numPredictors)), floor(numPredictors/3), floor(numPredictors/2), numPredictors))
numVarSplitVec <- c(floor(sqrt(numPredictors)):numPredictors)
## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)

## ############
## ## Uncomment the following to do cross-validation over and over.
## ## Number of times to do cross-validation.
## numCVs <- 100
## ## How many observations to reserve for testing each time.
## numLeaveOut <- round(0.10 * nrow(wideT))
## ## For matrix to hold cross-validation results.
## cvMSE <- matrix(NA, nrow(combos), ncol=numCVs)
## ############

## Number of times to do cross-validation.
numFolds <- 10
## Figure out which observations are in which fold by randomly
## assigning the numbers 1-10 to the various rows.  There are 93
## observations, so we assign an extra 1, 2, 3.
set.seed(425471)
numbersToAssign <- c(rep(1:10, 9), 1:3)
whichFold <- sample(numbersToAssign, replace=F)

## For matrix to hold cross-validation results.
cvMSE <- matrix(NA, nrow(combos), ncol=numFolds)


## Do cross-validation.
for (i in 1:numFolds){
  ## ############
  ## for (i in 1:numCVs){
  ## ##  Uncomment the following to do cross-validation over and over.
  ## ## Determine training and cross-validation set.
  ## whichLeaveOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)    
  ## subT <- wideT[-whichLeaveOut,] %>% select(-subj, -Rare, -days)
  ## cvsetT <- wideT[whichLeaveOut,] %>% select(-subj, -Rare, -days)
  ## ############

  ## Leave out all the rows assigned to i.
  whichLeaveOut <- whichFold==i
  cvsetT <- wideT[whichLeaveOut,] %>% select(-subj, -Rare, -days)
  subT <- wideT[!whichLeaveOut,] %>% select(-subj, -Rare, -days)
  
  for (j in 1:nrow(combos)){    
    rf <- randomForest(degdays ~ . , data=subT, mtry=combos[j,"numVarSplit"], ntree=combos[j,"numBtSamps"], importance=T)
    fitTest <- predict(rf, newdata=cvsetT)
    cvMSE[j,i] <- mean((fitTest - cvsetT$degdays)^2)
  }
}
rm(i, j)

combos$avgcvMSE <- apply(cvMSE, 1, mean)
## combos$q1cvMSE <- apply(cvMSE, 1, quantile, 0.25)
## combos$q3cvMSE <- apply(cvMSE, 1, quantile, 0.75)
ggplot(data=combos, aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) +
  geom_line()

## plot(range(numBtSampsVec),
##      range(as.vector(combos[,c("q1cvMSE", "q3cvMSE")])), type="n")
## my.colors <- c("black", "red", "green", "blue", "orange", "cyan", "magenta", "yellow", "brown")
## for (i in c(1:length(numVarSplitVec))){
##   subDF <- subset(combos, numVarSplit==numVarSplitVec[i])
##   with(subDF, lines(numBtSamps, avgcvMSE, col=my.colors[i]))
##   with(subDF, lines(numBtSamps, q1cvMSE, lty=2, col=my.colors[i]))
##   with(subDF, lines(numBtSamps, q3cvMSE, lty=2, col=my.colors[i]))  
## }
## rm(subDF, i)

yrng <- range(combos[,"avgcvMSE"])
ggplot(data=subset(combos, numVarSplit<=5), aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line() + ylim(yrng)
X11()
ggplot(data=subset(combos, (numVarSplit>5) & (numVarSplit<=8)), aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line() + ylim(yrng)
X11()
ggplot(data=subset(combos, numVarSplit>8), aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line() + ylim(yrng)
## mtry=6-8 and ntrees >= 3000
## ##################################################




## ##################################################
## Try random forests for classification.

## How many predictors are there?  (All taxa except "rare".)
numPredictors <- nrow(taxaT %>% filter(taxa!="Rare") %>% distinct(taxa))

## Try different numbers of bootstrap samples.
numBtSampsVec <- seq(500, 4000, by=500)
## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
## numVarSplitVec <- 6:8
## numVarSplitVec <- unique(c(floor(sqrt(numPredictors)), floor(numPredictors/3), floor(numPredictors/2), numPredictors))
numVarSplitVec <- c(floor(sqrt(numPredictors)):numPredictors)
## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)


## Number of times to do cross-validation.
numFolds <- 10
## Figure out which observations are in which fold by randomly
## assigning the numbers 1-10 to the various rows.  There are 93
## observations, so we assign an extra 1, 2, 3.
set.seed(246456)
numbersToAssign <- c(rep(1:10, 9), 1:3)
whichFold <- sample(numbersToAssign, replace=F)

## For matrix to hold fraction correctly classified.
errorFrac <- matrix(NA, nrow(combos), ncol=numFolds)


## Make a new version of the data frame with "days" variable as a
## factor, not a number.
factordaysT <- wideT
factordaysT$days <- as.factor(wideT$days)

## Do cross-validation.
for (i in 1:numFolds){

  ## Leave out all the rows assigned to i.
  whichLeaveOut <- whichFold==i
  cvsetT <- factordaysT[whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
  subT <- factordaysT[!whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
  
  for (j in 1:nrow(combos)){    
    rf <- randomForest(days ~ . , data=subT, mtry=combos[j,"numVarSplit"], ntree=combos[j,"numBtSamps"], importance=T)
    fitTest <- predict(rf, newdata=cvsetT)
    errorFrac[j,i] <- sum(fitTest != cvsetT$days)/length(fitTest)
  }
}
rm(i, j)

combos$avgerrorFrac <- apply(errorFrac, 1, mean)
ggplot(data=combos, aes(x=numBtSamps, y=avgerrorFrac, color=as.factor(numVarSplit))) +
  geom_line()


yrng <- range(combos[,"avgerrorFrac"])
ggplot(data=subset(combos, numVarSplit<=5), aes(x=numBtSamps, y=avgerrorFrac, color=as.factor(numVarSplit))) + geom_line() + ylim(yrng)
X11()
ggplot(data=subset(combos, (numVarSplit>5) & (numVarSplit<=8)), aes(x=numBtSamps, y=avgerrorFrac, color=as.factor(numVarSplit))) + geom_line() + ylim(yrng)
X11()
ggplot(data=subset(combos, numVarSplit>8), aes(x=numBtSamps, y=avgerrorFrac, color=as.factor(numVarSplit))) + geom_line() + ylim(yrng)
## mtry=4 and ntrees >= 2000
## ##################################################
