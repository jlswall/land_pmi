library("readxl")
library("tidyr")
library("tibble")
library("dplyr")
library("stringr")
library("ggplot2")
library("vegan")
library("colorspace")
library("randomForest")


## library("openxlsx")
## phylumAllT <- read.xlsx(fileNm, sheet="Phylum_all", skipEmptyCols=FALSE)
fileNm <- "../orig_data/Shane_all_skin_samples_taxo_bs_05_05_2017.xlsx"
rawAllT <- read_excel(path=fileNm, sheet="Phylum_all")

## The taxon column contains names of the form "p__name" or "k_name",
## and I will want to make a new column without the "p__" suffix
## later.  So, I rename this column now to avoid confusion later.
colnames(rawAllT)[colnames(rawAllT)=="taxon"] <- "taxonName"

## Some of these column names are duplicated.
## duplicated(colnames(phylumAllDF))
## This first happens at column 116.
## min(which(duplicated(colnames(phylumAllDF))))



## ##################################################
## Extract raw data for each pig in wide format.

## For now, I'm going to completely ignore columns 2 ("rankID") and 4
## ("daughterlevels").  I'm also going to ignore columsn beyond column
## "DJ" (column 114).  These appear to be summary columns.
## The columns I want to focus on are:
## column 3: "taxonName"
## column 5: "total" - to check that the other numbers are correct.
## columns 6-114
##    columns whose names start with "A#" - These are the observed
## counts for the individual pigs at the various times.
##    columns whose names start with "T#" - These are the averages
## of the observed counts over all pigs at the various times.  The
## names of these columns are in form T#days_#accumulatedDegreeDays.

## We're ignoring columns past 114.
mainT <- rawAllT[,1:114]


## Identify column names starting with "A".
namesA <- colnames(mainT)[substring(first=1, last=1, colnames(mainT))=="A"]
wideIndivT <- mainT[,c("taxonName", namesA)]


## Identify column names starting with "T".  The values in these
## columsn are the averages of the individual pigs at each time point.
## The names of these columns contain the number of days since death
## and the accumulated degree days.  The number of days since death
## immediately follows the "T", and the number of accumulated degree
## days follows the "_".
namesT <- colnames(mainT)[substring(first=1, last=1, colnames(mainT))=="T"]
## Separate the days from the accumulated degree days.
timeDF <- separate(data.frame(x=substring(namesT, first=2),
                              stringsAsFactors=F),
                   x, sep="_", into=c("days", "degdays"), convert=T)
## Note: number of days and accum. degree days are strongly correlated.
with(timeDF, cor(degdays, days))



## Extract the columns with the phylum names and the average counts
## across pigs for each time point.
wideAvgsT <- mainT[,c("taxonName", namesT)]

rm(namesA, namesT, mainT)
## ##################################################




## ##################################################
## Try to go from wide format to long format.

indivT <- wideIndivT %>%
  gather(indiv_time, counts, -taxonName) %>%
  separate(indiv_time, sep="_T", into=c("subj", "days"), convert=T)
  ## To add rows of missing values for combinations of days and
  ## subjects on which no samples are available (some subjects
  ## were not observed on certain days), add %>% and then this:
  ## complete(taxonName, days, subj)


## Next thing to check is whether we can get back the values in
## wideAvgsT when we do a summary for indivT.
chkAvgsT <- indivT %>%
  select(taxonName, days, counts) %>%
  group_by(taxonName, days) %>%
  summarize(avgs=mean(counts, na.rm=TRUE)) %>%
  spread(key=days,  value=avgs)
## Match the names to the timeDF frame.
matchNamesV <- na.omit(match(colnames(chkAvgsT), as.character(timeDF$days)))
chkNamesV <- paste(timeDF[matchNamesV,1], timeDF[matchNamesV,2], sep="_")
colnames(chkAvgsT) <- c("taxonName", paste0("T", chkNamesV))
## Order this in the same order as what I read in from the sheet.
reorderChkT <- chkAvgsT[match(wideAvgsT$taxonName, chkAvgsT$taxonName), colnames(wideAvgsT)]

## Compare with the averages I read in from the sheet.
all.equal(wideAvgsT[,1], reorderChkT[,1])  ## Ensure taxonNames in same order.
## Now check the counts, not the taxa names.
apply(wideAvgsT[,-1] - reorderChkT[,-1], 2, summary)
## There is a problem with column T1_27.  As an example, consider the
## row for p__Firmicutes.
subset(indivT, (days==1) & (taxonName=="p__Firmicutes"), "counts")
apply(subset(indivT, (days==1) & (taxonName=="p__Firmicutes"), "counts"), 2, mean)
subset(wideAvgsT, (taxonName=="p__Firmicutes"), "T1_27")
## It looks like they took the sum, not the mean.
apply(subset(indivT, (days==1) & (taxonName=="p__Firmicutes"), "counts"), 2, sum)
subset(wideAvgsT, (taxonName=="p__Firmicutes"), "T1_27")
## Look at all the values for "T1_27".
## cbind(reorderChkT[,"T1_27"], wideAvgsT[,"T1_27"], reorderChkT[,"T1_27"]- wideAvgsT[,"T1_27"] )

rm(chkAvgsT, matchNamesV, chkNamesV)
## ##################################################



## ##################################################
## Make other adjustments to the dataset so that it's easier to use.

## Check that all the phylum counts add up to the k__Bacteria counts.
tmp1 <- indivT %>%
  filter(taxonName!="k__Bacteria") %>%
  group_by(days, subj) %>%
  summarize(counts=sum(counts))
tmp2 <- subset(indivT, taxonName=="k__Bacteria", c("days", "subj", "counts"))
if (!all.equal(tmp1, tmp2))
  stop("Phylum counts don't add up to k__Bacteria counts")
rm(tmp1, tmp2)


## Remove the k__Backteria row, since it is just the totals of the
## phyla.  Also, include accum. degree days in the tibble.
indivT <- indivT %>%
  filter(taxonName!="k__Bacteria") %>%
  left_join(timeDF, by="days")


## Make a new, more readable taxa column by stripping the "p__"
## from the beginning of each taxon name.
indivT$taxa <- gsub(indivT$taxonName, pattern="p__", replacement="")
## The column name "[Thermi]" causes problems for functions expecting
## traditional data frame column names.
indivT$taxa <- gsub(indivT$taxa, pattern="\\[", replacement="")
indivT$taxa <- gsub(indivT$taxa, pattern="]", replacement="")
## Remova the taxonName column from the tibble to avoid confusion with
## the next taxa column.
indivT <- indivT %>% select(-taxonName)
## ##################################################



## ##### WORKING HERE: Would like to move this part down so that we'll
## be working with just the most common taxa.


## ##################################################
## For each individual pig and each day, find the percentage of all
## counts represented by each taxa.

ctsByDaySubjT <- indivT %>%
  group_by(days, degdays, subj) %>%
  summarize(totals = sum(counts))

indivT <- indivT %>%
  left_join(ctsByDaySubjT) %>%
  mutate(fracByDaySubj=counts/totals) %>%
  select(-totals)
## ##################################################




## ##################################################
## Some taxa don't occur frequently.  We identify these, classify them
## as rare, and then reformat the data accordingly.
## Pechal et al (2013) call "rare" any taxa with <3% of total
## abundance.

## Find taxa which make up >= 3% of the total counts, over all days
## and subjects.  These are the "not rare" taxa, which I'll call
## "common".  Don't include the "unclassified" category of taxa, as
## unclassified taxa counts/fractions won't be used in the individual
## models.
commonByTotalV <- unlist(indivT %>%
    group_by(taxa) %>%
    summarize(taxatotal = sum(counts)) %>%
    mutate(frac = taxatotal/sum(taxatotal)) %>%
    filter(frac >= 0.03) %>%
    filter(taxa != "unclassified") %>%
    select(taxa)
    )
## See a list of all taxa percentages sorted in descending order:
## data.frame( indivT %>%
##     group_by(taxa) %>%
##     summarize(taxatotal = sum(counts)) %>%
##     mutate(frac = taxatotal/sum(taxatotal)) %>%
##     arrange(desc(frac))
##     )


## Taxa counts can vary widely by day and individual.  So, another way
## to figure out which taxa are "common" would be to include any taxa
## making up at least 3% of the total count on at least one particular
## day for any particular subject.  As above, we don't include the
## "unclassified" category of taxa.  (For phylum taxa, unclassified
## does have a maximum fraction of about 0.0343 on day 47, subject A6.
## For other days and subjects, it's less than 0.02).
commonByDaySubjV <- unlist(indivT %>%
                           group_by(taxa) %>%
                           summarize(maxFracByDaySubj = max(fracByDaySubj)) %>%
                           filter(maxFracByDaySubj >= 0.03) %>%
                           filter(taxa != "unclassified") %>%
                           select(taxa)
                           )
## See a list of maximum taxa percentages sorted in descending order:
## indivT %>%
##   group_by(taxa) %>%
##   summarize(maxFracByDaySubj = max(fracByDaySubj)) %>%
##   arrange(desc(maxFracByDaySubj))


## Check whether the "common" taxa are the same whether we use the 3%
## cutoff as 3% of the taxa totals or 3% of the totals by day and
## subject.
{
if (identical(commonByTotalV, commonByDaySubjV))
  commonTaxaNamesV <- commonByTotalV
else
  stop("Taxa common enough to be included in further analyses are different when calculating for totals vs. for day/subj.")
}
rm(commonByTotalV, commonByDaySubjV)



## Rename unclassified taxa and taxa that occur less than 3% of the
## time to "rare or uncl.".  Then, sum the counts for all these
## rare or uncl. taxa into one row.
commontaxaT <- indivT
commontaxaT[!(commontaxaT$taxa %in% commonTaxaNamesV), "taxa"] <- "Rare or uncl."
commontaxaT <- commontaxaT %>%
  group_by(days, degdays, subj, taxa) %>%
  summarize(counts = sum(counts), fracByDaySubj=sum(fracByDaySubj))
rm(commonTaxaNamesV)
## ##################################################




## ##################################################
## Treat each pig as a "site"/"community" and each day as a time step.
## The repsonse is multivariate because there are counts for many
## different phyla.

## We need to get the data into "community matrix format", which is
## basically a wide format.
wideT <- commontaxaT %>%
##  wideT <- indivT %>%
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
## Graphics illustrating the variability in counts, by individual and
## by day.

## Counts are highly variable, by individual and by day.
ggplot(ctsByDaySubjT, aes(degdays, totals)) +
  geom_point(aes(color=subj)) +
  xlab("Degree days") +
  ylab("Total taxa counts by degree day and subject") +
  labs(color="Subject")

## For each bacteria, plot counts for each pig vs. accum. degree days.
## Using raw counts, we can see the large variability among individuals.
ggplot(commontaxaT, aes(degdays, counts)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa) +
  xlab("Degree days") +
  ylab("Counts by degree day and subject") +
  labs(color="Subject")


## Make the stacked bar chart using raw counts by degree days.
ggplot(commontaxaT, aes(degdays)) +
  geom_bar(aes(weight=counts, fill=taxa)) +
  xlab("Degree days") +
  ylab("Counts by degree day and taxa, combined over subjects")


## Assess variability among pigs on the first few days, in terms of counts.
ggplot(subset(commontaxaT, days <= 5), aes(x=subj, y=counts, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days)
## ##################################################



## ##################################################
## Bar charts using fractions associated with each taxa, rather than
## raw counts.

## Summarize the counts by day and by day-taxa, calculate percentages
## of each taxa per day (across all pigs).
bydayT <- commontaxaT %>%
  group_by(days, degdays, taxa) %>%
  summarize(ctsByTaxaDay = sum(counts)) %>%
  left_join(commontaxaT %>%
              group_by(degdays) %>%
              summarize(ctsByDay = sum(counts)),
            by = "degdays"
  ) %>%
  mutate(fracByDayTaxa = ctsByTaxaDay / ctsByDay)


## Make stacked bar chart using percentages of taxa for each degree
## day based on the total count for that day across all pigs.
ggplot(bydayT, aes(degdays)) +
  geom_bar(aes(weight=fracByDayTaxa, fill=taxa)) +
  labs(x="Accum. degree days", y="Relative abundance")

## Make stacked bar chart to compare with Figure 1 in Pechal et al
## (2013).
ggplot(subset(bydayT, days <= 5), aes(days)) +
  geom_bar(aes(weight=fracByDayTaxa, fill=taxa)) +
  labs(x="Days", y="Relative abundance")


## Assess variability among pigs on the first few days, using fractions.
ggplot(subset(commontaxaT, days <= 5),
       aes(x=subj, y=fracByDaySubj, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days)
## ##################################################




## ##################################################
## Try random forests.

## Move back to wide format.
widePercT <- commontaxaT %>%
  ungroup() %>%
  select(degdays, subj, taxa, fracByDaySubj) %>%
  spread(taxa, fracByDaySubj)
## One group ("rare or uncl.") has spaces in it, which can cause
## problems when the tibble is converted to a data.frame.
colnames(widePercT) <- gsub(colnames(widePercT), pattern=" ", replacement="_")
    

## Pick subset of the data to train on.
trainingIndices <- sort(sample(1:nrow(widePercT), size=74, replace=F))
rf <- randomForest(degdays ~ . -subj -Rare_or_uncl., data=widePercT[trainingIndices,], importance=T)
imp.rf <- importance(rf)
varImpPlot(rf)
## For phylum taxa: Firmicutes and Actinobacteria strongest.  Less
## strong are Proteobacteria and Bacteroidetes.

## Now try to predict for those observations not in the training set.
yhatTest <- predict(rf, newdata=widePercT[-trainingIndices,])
mean((yhatTest - widePercT[-trainingIndices,"degdays"])^2)

impvar <- rownames(imp.rf)[order(imp.rf[,1], decreasing=T)]
for (i in seq_along(impvar))
  partialPlot(rf, pred.data=as.data.frame(widePercT[trainingIndices,]), x.var=impvar[i], main=paste("Partial Dependence on", impvar[i]))
rm(i)
## ##################################################




## ##################################################
## Try out gam models.

## For each bacteria, plot counts for each pig vs. accum. degree days.
## Using raw counts, we can see the large variability among individuals.
ggplot(commontaxaT %>% filter(taxa!="Rare or uncl."),
      aes(degdays, fracByDaySubj)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa) +
  labs(x="Degree days", y="Fraction by degree day and subject") +
  labs(color="Subject")

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
