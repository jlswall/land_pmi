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
## Try random forests.

## How many predictors are there?  (All taxa except "rare".)
numPredictors <- nrow(taxaT %>% filter(taxa!="Rare") %>% distinct(taxa))

## Try different numbers of bootstrap samples.
numBtSamps <- seq(250, 3000, by=250)
## numBtSamps <- seq(500, 5000, by=500)
## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
numVarSplit <- unique(c(2, floor(sqrt(numPredictors)), floor(numPredictors/3), floor(numPredictors/2), numPredictors))
## numVarSplit <- c(floor(sqrt(numPredictors)):numPredictors)
## Form matrix with all combinations of these.
combos <- expand.grid(numBtSamps=numBtSamps, numVarSplit=numVarSplit)
## Add another column to contain the cross-validation MSE.
combos$cvMSE <- NA

## Number of times to do cross-validation.
numCVs <- 100

set.seed(347194)
cvMSE <- rep(NA, length(combos))
for (i in 1:nrow(combos)){
  
  ## ## Do cross-validation, leaving one pig out at a time.
  ## foldMSE <- NULL
  ## for (pig.k in unique(wideT$subj)){
    
  ##   subT <- wideT %>% filter(subj!=pig.k) %>% select(-subj, -Rare, -degdays)
  ##   cvsetT <- wideT %>% filter(subj==pig.k) %>% select(-subj, -Rare, -degdays)
  ##   rf <- randomForest(days ~ . , data=subT, mtry=combos[i,"numVarSplit"], ntree=combos[i,"numBtSamps"], importance=T)
  ##   fitTest <- predict(rf, newdata=cvsetT)
  ##   foldMSE <- c(foldMSE, sum((fitTest - cvsetT$days)^2)/length(fitTest))
  ## }
  
  ## Do cross-validation, leaving one about 10% of the data at time.
  foldMSE <- NULL
  numLeaveOut <- round(0.10 * nrow(wideT))
  for (j in 1:numCVs){
    
    whichLeaveOut <- sample(1:nrow(wideT), size=numLeaveOut, replace=F)    
    subT <- wideT[-whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
    cvsetT <- wideT[whichLeaveOut,] %>% select(-subj, -Rare, -degdays)
    
    rf <- randomForest(days ~ . , data=subT, mtry=combos[i,"numVarSplit"], ntree=combos[i,"numBtSamps"], importance=T)
    fitTest <- predict(rf, newdata=cvsetT)
    foldMSE <- c(foldMSE, mean((fitTest - cvsetT$days)^2))
  }
  
  combos$cvMSE[i]<- mean(foldMSE)
}
rm(i, foldMSE)

ggplot(data=combos, aes(x=numBtSamps, y=sqrt(cvMSE), color=as.factor(numVarSplit))) +
  geom_line()
## ##################################################



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
