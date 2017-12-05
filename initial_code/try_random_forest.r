library("tidyverse")
library("randomForest")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "orders"

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
## Try random forests for regression.

## How many predictors are there?  (All taxa except "rare".)
numPredictors <- nrow(taxaT %>% filter(taxa!="Rare") %>% distinct(taxa))

## Try different numbers of bootstrap samples.
numBtSampsVec <- seq(500, 4000, by=500)
## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
numVarSplitVec <- 6:8
## numVarSplitVec <- unique(c(floor(sqrt(numPredictors)), floor(numPredictors/3), floor(numPredictors/2), numPredictors))
## numVarSplitVec <- c(floor(sqrt(numPredictors)):numPredictors)
## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)

## Number of times to do cross-validation.
numCVs <- 50
## For matrix to hold cross-validation results.
cvMSE <- matrix(NA, nrow(combos), ncol=numCVs)


set.seed(347194)
## Do cross-validation, leaving one about 10% of the data at time.
numLeaveOut <- round(0.10 * nrow(wideT))
for (i in 1:numCVs){

  ## Determine training and cross-validation set.
  whichLeaveOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)    
  subT <- wideT[-whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
  cvsetT <- wideT[whichLeaveOut,] %>% select(-subj, -Rare, -degdays)

  for (j in 1:nrow(combos)){    
    rf <- randomForest(days ~ . , data=subT, mtry=combos[j,"numVarSplit"], ntree=combos[j,"numBtSamps"], importance=T)
    fitTest <- predict(rf, newdata=cvsetT)
    cvMSE[j,i] <- mean((fitTest - cvsetT$days)^2)
  }
}
rm(i, j)

avgcvMSE <- apply(cvMSE, 1, mean)
combos$avgcvMSE <- apply(cvMSE, 1, mean)
combos$q1cvMSE <- apply(cvMSE, 1, quantile, 0.25)
combos$q3cvMSE <- apply(cvMSE, 1, quantile, 0.75)
ggplot(data=combos, aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) +
  geom_line()

plot(range(numBtSampsVec),
     range(as.vector(combos[,c("q1cvMSE", "q3cvMSE")])), type="n")
my.colors <- c("black", "red", "green", "blue", "orange", "cyan", "magenta", "yellow", "brown")
for (i in c(1:length(numVarSplitVec))){
  subDF <- subset(combos, numVarSplit==numVarSplitVec[i])
  with(subDF, lines(numBtSamps, avgcvMSE, col=my.colors[i]))
  with(subDF, lines(numBtSamps, q1cvMSE, lty=2, col=my.colors[i]))
  with(subDF, lines(numBtSamps, q3cvMSE, lty=2, col=my.colors[i]))  
}
rm(subDF, i)

ggplot(data=subset(combos, (numVarSplit<10) & (numVarSplit>4)), aes(x=numBtSamps, y=avgcvMSE, color=as.factor(numVarSplit))) + geom_line()
## ##################################################



## ##################################################
## Try random forests for classification.

## How many predictors are there?  (All taxa except "rare".)
numPredictors <- nrow(taxaT %>% filter(taxa!="Rare") %>% distinct(taxa))

## Try different numbers of bootstrap samples.
numBtSampsVec <- seq(50, 500, by=50)
## numBtSamps <- seq(500, 5000, by=500)
## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
numVarSplitVec <- unique(c(floor(sqrt(numPredictors)), floor(numPredictors/3), floor(numPredictors/2), numPredictors))
## numVarSplit <- c(floor(sqrt(numPredictors)):numPredictors)
## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSampsVec, numVarSplit=numVarSplitVec)
## Add another column to contain the cross-validation MSE.
combos$avgcvMSE <- NA
combos$q1cvMSE <- NA
combos$q3cvMSE <- NA

## Number of times to do cross-validation.
numCVs <- 50


set.seed(431890)
cvMSE <- rep(NA, length(combos))
for (i in 1:nrow(combos)){
  
  ## Do cross-validation, leaving one about 10% of the data at time.
  foldMSE <- NULL
  numLeaveOut <- round(0.10 * nrow(wideT))
  for (j in 1:numCVs){
    
    whichLeaveOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)    
    subT <- wideT[-whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
    cvsetT <- wideT[whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
    
    rf <- randomForest(as.factor(days) ~ . , data=subT, mtry=combos[i,"numVarSplit"], ntree=combos[i,"numBtSamps"], importance=T)
    fitTest <- predict(rf, newdata=cvsetT)
    foldMSE <- c(foldMSE, mean((fitTest - cvsetT$days)^2))
  }
  
  combos$avgcvMSE[i]<- mean(foldMSE)
  combos$q1cvMSE[i]<- quantile(foldMSE, 0.25)
  combos$q3cvMSE[i]<- quantile(foldMSE, 0.75)
}
rm(i, foldMSE)






## ##################################################
## Try random forests.

## Pick subset of the data to train on.
trainingIndices <- sort(sample(1:nrow(wideT), size=60, replace=F))
rf <- randomForest(degdays ~ . -subj -Rare, data=wideT[trainingIndices,], importance=T)
imp.rf <- importance(rf)
varImpPlot(rf)
## For phylum taxa: Firmicutes strongest.  Bacteroidetes not strong.

## Now try to predict for those observations not in the training set.
yhatTest <- predict(rf, newdata=wideT[-trainingIndices,])
mean((yhatTest - wideT[-trainingIndices,"degdays"])^2)

par(mfrow=c(3,4))
impvar <- rownames(imp.rf)[order(imp.rf[,1], decreasing=T)]
for (i in seq_along(impvar))
  partialPlot(rf, pred.data=as.data.frame(wideT[trainingIndices,]), x.var=impvar[i], main=paste("Partial Dependence on", impvar[i]))
rm(i)
## ##################################################



## ##################################################
## Try out gam models.


## Try gamsel package.
trainX <- as.data.frame(wideT %>%
                        select(-one_of("degdays", "subj", "Rare_or_uncl."))
                        )[trainingIndices,]
trainY <- as.matrix(as.data.frame(wideT[trainingIndices, "degdays"]))
gamsel.out <- gamsel(x=trainX, y=trainY, dfs=rep(10, ncol(trainX)))
summary(gamsel.out)
gamsel.cv <- cv.gamsel(x=trainX, y=trainY, dfs=rep(10, ncol(trainX)))


gam.fit <- gam(degdays ~ s(Firmicutes) + s(Actinobacteria), data=wideT[trainingIndices,])
big.fit <- gam(degdays ~ s(Firmicutes) + s(Actinobacteria) + s(Proteobacteria), data=wideT[trainingIndices,])
summary(gam.fit)
yhatTest <- predict(gam.fit, newdata=wideT[-trainingIndices,])

## Show fitted values vs. real values.
plot(unlist(wideT[trainingIndices, "degdays"]), gam.fit$fitted)
cor(unlist(wideT[trainingIndices, "degdays"]), gam.fit$fitted)
## Show predicted values vs. real values.
plot(unlist(wideT[-trainingIndices, "degdays"]), unlist(yhatTest))
cor(unlist(wideT[-trainingIndices, "degdays"]), unlist(yhatTest)
## ##################################################



