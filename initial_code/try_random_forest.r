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


## Number of possible predictors.
numPredictors <- nrow(taxaT %>% filter(taxa!="Rare") %>% distinct(taxa))

## Try different numbers of bootstrap samples.
numBtSamps <- seq(500, 5000, by=500)

## Try different values for mtry (which represents how many variables
## can be chosen from at each split of the tree).
numVarSplit <- c(floor(sqrt(numPredictors)):numPredictors)
set.seed(347194)
cvMSE <- matrix(NA, nrow=length(numBtSamps), ncol=length(numVarSplit))
for (i in 1:length(numBtSamps)){

  ntrees <- numBtSamps[i]
  
  for (j in 1:length(numVarSplit)){

    m <- numVarSplit[j]
  
    ## Do cross-validation, leaving one pig out at a time.
    cvMSEs <- NULL
    for (pig.k in unique(wideT$subj)){
      
      subT <- wideT %>% filter(subj != pig.k) %>% select(-subj, -Rare, -degdays)
      cvsetT <- wideT %>% filter(subj == pig.k) %>% select(-subj, -Rare, -degdays)
      rf <- randomForest(days ~ . , data=subT, mtry=m, ntree=ntrees, importance=T)
      cvTest <- predict(rf, newdata=cvsetT)
      cvMSEs <- c(cvMSEs, sum((cvTest - cvsetT$days)^2)/length(cvTest))
    }
    cvMSE[i, j]<- mean(cvMSEs)
  }
}
rm(i, pig.k, j)
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
