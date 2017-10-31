library("tidyverse")


## Read in cleaned-up phyla taxa.
phylumT <- read_csv("phyla_massaged.csv")
## Read in cleaned-up families taxa.
familyT <- read_csv("families_massaged.csv")



## ##################################################
## Move to wide format.
wideT <- phylumT %>%
  filter(taxa!="Rare") %>%
  select(days, degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay) 
## ##################################################


## ##################################################
## At random, select observations to make up the training dataset.

set.seed(345948)
trainingIndices <- sort(sample(1:nrow(wideT), size=60, replace=F))
trainT <- data.frame(wideT[trainingIndices,])
rownames(trainT) <- paste0("D", trainT$days, "_S", trainT$subj )
## ##################################################


## ##################################################
## Do PCA on training dataset.
## References for PCA in R:
## https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/
## Pages 402-403 of "An Introduction to Statistical Learning"

pr.out <- prcomp(trainT %>% select(-days, -degdays, -subj), scale=T)
biplot(pr.out, scale=0, cex=0.7)

## Fraction of variance explained
fracVarExp <- pr.out$sdev^2/sum(pr.out$sdev^2)

## Bar chart of variance explained.
pdf(file="phyla_prcomp_expl_barchart.pdf")
barplot(cumsum(fracVarExp), names.arg=paste0("PC ", 1:length(pr.out$sdev)), ylab="Fraction of variance explained")
abline(h=0.75, lty=2)
dev.off()
## ##################################################




## Pick subset of the data to train on.
trainingIndices <- sort(sample(1:nrow(wideT), size=60, replace=F))
rf <- randomForest(degdays ~ . -subj -Rare, data=wideT[trainingIndices,], importance=T)
imp.rf <- importance(rf)
varImpPlot(rf)
## For phylum taxa: Firmicutes strongest.  Bacteroidetes not strong.

## Now try to predict for those observations not in the training set.
yhatTest <- predict(rf, newdata=wideT[-trainingIndices,])
mean((yhatTest - wideT[-trainingIndices,"degdays"])^2)

par(mfrow=c(2,2))
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
