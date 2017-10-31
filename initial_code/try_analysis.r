library("tidyverse")
library("vegan")
library("randomForest")
library("colorspace")


## Read in cleaned-up phyla taxa.
phylumT <- read_csv("phyla_massaged.csv")
## Read in cleaned-up families taxa.
familyT <- read_csv("families_massaged.csv")


## ##################################################
## Treat each pig as a "site"/"community" and each day as a time step.
## The repsonse is multivariate because there are counts for many
## different phyla.

## We need to get the data into "community matrix format", which is
## basically a wide format.
wideT <- phylumT %>%
  select(days, degdays, subj, taxa, counts) %>%
  spread(taxa, counts)
## Now, break this into 2 pieces.  One piece is the "community
## matrix", which has a column for each phylum and with the count for
## that phylum for a given day and subject.  The second piece is the
## covariates (subject and day).
communityCovariates <- wideT[,1:3]
communityCounts <- wideT[,-c(1:3)]


## Fit the PERMANOVA model
## First without strata:
## trynostrata <- adonis(communityCounts ~ as.factor(days), data=communityCovariates, permutations=5000)
## Then, with strata.  (See
## http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html).
trystrata <- adonis(communityCounts ~ as.factor(days), data=communityCovariates, strata=communityCovariates$subj, permutations=5000)


## Try using nonmetric multidimensional scaling.
trymds <- metaMDS(communityCounts, k=2, trymax=1000)
stressplot(trymds)
plot(trymds)
ordiplot(trymds,type="n")
orditorp(trymds, display="species", col="red", air=0.01,
         labels=make.cepnames(names(communityCounts)))
## Make colors.
myColors <- diverge_hcl(length(unique(communityCovariates$days)))
orditorp(trymds, display="sites", air=0.01, col=myColors[as.numeric(as.factor(communityCovariates$days))], labels=as.character(communityCovariates$days))
## ordihull(trymds, groups=communityCovariates$days, draw="polygon", col=myColors)
rm(myColors)

rm(trystrata, trymds)
## ##################################################



## ##################################################
## Scatterplot of percentage composition by degree day, with each
## bacteria a different color.  No distinguishing between individuals.

## For phylum taxa:
ggplot(phylumT %>% filter(taxa!="Rare"),
      aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=taxa)) +
  labs(x="Degree days", y="Composition fraction") +
  labs(color="Subject")
ggsave("phyla_scatter_frac_by_degday.pdf", width=7, height=4, units="in")

## For family taxa:
ggplot(familyT %>% filter(taxa!="Rare"),
      aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=taxa)) +
  labs(x="Degree days", y="Composition fraction") +
  labs(color="Subject")
ggsave("families_scatter_frac_by_degday.pdf", width=7, height=4, units="in")
## ##################################################



## ##################################################
## For each bacteria, scatterplot of fractions vs. accum. degree days
## (for each pig).

## For phylum taxa:
ggplot(phylumT %>% filter(taxa!="Rare"),
      aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa) +
  labs(x="Degree days", y="Composition fraction") +
  labs(color="Subject")

## For family taxa:
ggplot(familyT %>% filter(taxa!="Rare"),
      aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa) +
  labs(x="Degree days", y="Composition fraction") +
  labs(color="Subject")
ggsave("families_scatter_frac_by_degday_by_taxa.pdf", width=7, height=, units="in")
## ##################################################



## ##################################################
## Try random forests.

## Move back to wide format.
widePercT <- familyT %>%
  ungroup() %>%
  select(degdays, subj, taxa, fracBySubjDay) %>%
  spread(taxa, fracBySubjDay)
## One group ("rare or uncl.") has spaces in it, which can cause
## problems when the tibble is converted to a data.frame.
## colnames(widePercT) <- gsub(colnames(widePercT), pattern=" ", replacement="_")
    

## Pick subset of the data to train on.
trainingIndices <- sort(sample(1:nrow(widePercT), size=60, replace=F))
rf <- randomForest(degdays ~ . -subj -Rare, data=widePercT[trainingIndices,], importance=T)
imp.rf <- importance(rf)
varImpPlot(rf)
## For phylum taxa: Firmicutes strongest.  Bacteroidetes not strong.

## Now try to predict for those observations not in the training set.
yhatTest <- predict(rf, newdata=widePercT[-trainingIndices,])
mean((yhatTest - widePercT[-trainingIndices,"degdays"])^2)

par(mfrow=c(2,2))
impvar <- rownames(imp.rf)[order(imp.rf[,1], decreasing=T)]
for (i in seq_along(impvar))
  partialPlot(rf, pred.data=as.data.frame(widePercT[trainingIndices,]), x.var=impvar[i], main=paste("Partial Dependence on", impvar[i]))
rm(i)
## ##################################################



## ##################################################
## Try out gam models.


## Try gamsel package.
trainX <- as.data.frame(widePercT %>%
                        select(-one_of("degdays", "subj", "Rare_or_uncl."))
                        )[trainingIndices,]
trainY <- as.matrix(as.data.frame(widePercT[trainingIndices, "degdays"]))
gamsel.out <- gamsel(x=trainX, y=trainY, dfs=rep(10, ncol(trainX)))
summary(gamsel.out)
gamsel.cv <- cv.gamsel(x=trainX, y=trainY, dfs=rep(10, ncol(trainX)))


gam.fit <- gam(degdays ~ s(Firmicutes) + s(Actinobacteria), data=widePercT[trainingIndices,])
big.fit <- gam(degdays ~ s(Firmicutes) + s(Actinobacteria) + s(Proteobacteria), data=widePercT[trainingIndices,])
summary(gam.fit)
yhatTest <- predict(gam.fit, newdata=widePercT[-trainingIndices,])

## Show fitted values vs. real values.
plot(unlist(widePercT[trainingIndices, "degdays"]), gam.fit$fitted)
cor(unlist(widePercT[trainingIndices, "degdays"]), gam.fit$fitted)
## Show predicted values vs. real values.
plot(unlist(widePercT[-trainingIndices, "degdays"]), unlist(yhatTest))
cor(unlist(widePercT[-trainingIndices, "degdays"]), unlist(yhatTest)
## ##################################################
