library("tidyverse")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "orders"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0(taxalevel, "_massaged.csv"))
## ##################################################


## ##################################################
## Move to wide format.
wideT <- taxaT %>%
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
pdf(file=paste0(taxalevel, "_prcomp_expl_barchart.pdf"))
barplot(cumsum(fracVarExp), names.arg=paste0("PC ", 1:length(pr.out$sdev)), ylab="Fraction of variance explained")
abline(h=0.75, lty=2)
dev.off()
## ##################################################
