library("readxl")
library("tidyr")
library("tibble")
library("dplyr")
library("stringr")
library("ggplot2")
library("randomForest")


## library("openxlsx")
## phylumAllT <- read.xlsx(fileNm, sheet="Phylum_all", skipEmptyCols=FALSE)
fileNm <- "../orig_data/Shane_all_skin_samples_taxo_bs_05_05_2017.xlsx"
rawAllT <- read_excel(path=fileNm, sheet="family_all_bsedit")

## The "Family" column contains some names that are strange because
## they are surrounded by square brackets (e.g., [Tissierellaceae]").
## I will want to make a new column later with only text-based names.
## So, I rename this column now to avoid confusion later.
colnames(rawAllT)[colnames(rawAllT)=="Family"] <- "origName"

## Column names are unique in this sheet (unlike the "Phylum_all"
## sheet.)
sum(duplicated(colnames(rawAllT)))



## ##################################################
## Concentrate on the counts for the individual pigs on the various
## days.  Transfer from wide format to long format.  Check that the
## columns whose names start with "T" are actually the averages I
## think they are, and I'll check that totals column and row seems
## correct.

## Column 1: taxa name (last row, labeled "Bacteria", contains
##    totals over all taxa)
## Columns 2-110: 
##    columns whose names start with "A" - These are the observed
## counts for the individual pigs at the various times.
##    columns whose names start with "T" - These are the averages
## of the observed counts over all pigs at the various times.  The
## names of these columns are in form T#days_#accumulatedDegreeDays.
## Column 111: "total" - Want to check that this is correct.
## Columns >111:  seem to be summary columns, will check later.

## For this part, we're working with just the first 111 columns.
mainT <- rawAllT[,1:111]


## #######################
## Put individual counts and average counts into different tables.

## Identify column names starting with "A". Save these as the counts
## for individual pigs on the various days.
namesA <- colnames(mainT)[substring(first=1, last=1, colnames(mainT))=="A"]
wideIndivT <- mainT[,c("origName", namesA)]

## Identify column names starting with "T".  The values in these
## columns are the averages of the individual pigs at each time point.
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


## Extract the columns with the taxa names and the average counts
## across pigs for each time point.
wideAvgsT <- mainT[,c("origName", namesT)]

rm(namesA, namesT)
## #######################


## #######################
## Go from wide format to long format.  Check the averages and totals
## columns.

## Go from wide to long format.
indivT <- wideIndivT %>%
  gather(indiv_time, counts, -origName) %>%
  separate(indiv_time, sep="_T", into=c("subj", "days"), convert=T)
  ## To add rows of missing values for combinations of days and
  ## subjects on which no samples are available (some subjects
  ## were not observed on certain days), add %>% and then this:
  ## complete(taxonName, days, subj)


## Check whether we can get back the values in wideAvgsT when we do a
## summary for indivT.
chkAvgsT <- indivT %>%
  select(origName, days, counts) %>%
  group_by(origName, days) %>%
  summarize(avgs=mean(counts)) %>%
  spread(key=days,  value=avgs)
## Match the names to the timeDF frame.
matchNamesV <- na.omit(match(colnames(chkAvgsT), as.character(timeDF$days)))
chkNamesV <- paste(timeDF[matchNamesV,1], timeDF[matchNamesV,2], sep="_")
colnames(chkAvgsT) <- c("origName", paste0("T", chkNamesV))
## Order this in the same order as what I read in from the sheet.
reorderChkT <- chkAvgsT[match(wideAvgsT$origName, chkAvgsT$origName), colnames(wideAvgsT)]


## Compare with the averages I read in from the sheet.
## First, ensure taxonNames in same order.
all.equal(wideAvgsT[,1], reorderChkT[,1])
## Now check the counts, not the taxa names.
apply(wideAvgsT[,-1] - reorderChkT[,-1], 2, summary)
## There is a problem with column T1_27.  As an example, consider the
## counts for Clostridiaceae on that day.
subset(indivT, (days==1) & (origName=="Clostridiaceae"), "counts")
## The average is given by:
apply(subset(indivT, (days==1) & (origName=="Clostridiaceae"), "counts"), 2, mean)
## But, this is not the average we read from the sheet.
subset(wideAvgsT, (origName=="Clostridiaceae"), "T1_27")
## It looks like they took the sum, not the mean.
apply(subset(indivT, (days==1) & (origName=="Clostridiaceae"), "counts"), 2, sum)
## Look at all the values for "T1_27".
## cbind(reorderChkT[,"T1_27"], wideAvgsT[,"T1_27"], reorderChkT[,"T1_27"]- wideAvgsT[,"T1_27"] )
rm(chkAvgsT, matchNamesV, chkNamesV)
## #######################


## #######################
## Check that the total counts for each taxa match the "total" column
## (column #111).  There are a lot of these to check, so we take the
## absolute differences between the our total counts and the "total"
## column and make sure that the biggest difference is 0.
max(indivT %>%
  group_by(origName) %>%
  summarize(totalCt = sum(counts)) %>%
  left_join(rawAllT %>% select(origName, total)) %>%
  mutate(absDiffOrigMyCalc = abs(total - totalCt)) %>%
  select(absDiffOrigMyCalc)
  )
## #######################


## #######################
## Check that the total counts for each subject on each day match the
## "Bacteria" row (row #233 in the tibble, #234 in the worksheet).

## Save these totals in a table for use in calculating percentages.
totalCtbySubjDay <- indivT %>%
  filter(origName!="Bacteria") %>%
  group_by(days, subj) %>%
  summarize(total=sum(counts))

## Compare the totals calculated above with the last row ("Bacteria")
## of the individual counts.
all.equal(
    unite(totalCtbySubjDay, subj_day, subj, days, sep="_T") %>%
      spread(key=subj_day, value=total),
    wideIndivT %>% filter(origName=="Bacteria") %>% select(-origName)
)
## #######################

rm(mainT, wideAvgsT, wideIndivT, reorderChkT)
## ##################################################




## ##################################################
## Columns 113-223 ("DI"-"HO") appear to be percentages for each taxa, by
## day and pig.  Check these.
## Columns 225-242: I haven't checked these yet.

## #######################
## Organize the information.

## Put these columns in their own table.
widePercT <- rawAllT[,113:223]

## The first column contains the family names.
colnames(widePercT)[1] <- "origName"

## Identify column names starting with "A" (individuals A1-A6).
namesA <- colnames(widePercT)[substring(first=1, last=1, colnames(widePercT))=="A"]
wideIndivPercT <- rawAllT[,c("origName", namesA)]


## Identify column names starting with "T" (averages across
## individuals).
namesT <- colnames(widePercT)[substring(first=1, last=1, colnames(widePercT))=="T"]
wideAvgsPercT <- rawAllT[,c("origName", namesA)]
## #######################


## #######################
## Try to take the individual percentages from wide format to long
## format.
## indivPercT <- wideIndivPercT %>%
##   gather(indiv_time, perc, -origName) %>%
##   separate(indiv_time, sep="_T", into=c("subj", "days_with_extra"), convert=T) %>%
##   separate(days_with_extra, sep="__", into=c("days", "extra_stuff"), convert=T) %>%
##   select(-extra_stuff)
## #######################


## #######################
## Check that these percentages straight from the worksheet are the
## same as what we calculate based on the individual counts.

## First, I calculate these percentages based on the individuals
## counts and the sums (based on those counts) that I calculated
## earlier.  I add this column to the main table.
indivT <- indivT %>%
  left_join(totalCtbySubjDay) %>%
  mutate(percByDaySubj = 100*counts/total) %>%
  select(-total)


## Try to put these percentages in wide format for comparison with the
## raw numbers from the worksheet.
chkPercT <- indivT %>%
  select(-counts) %>%
  mutate(extrachar="1") %>%
  unite(subj_days, subj, days, sep="_T") %>%
  unite(subj_days, subj_days, extrachar, sep="__") %>%
  spread(key=subj_days, value=percByDaySubj)

## Now check to see if this matches the numbers we got from the spreadsheet.
if (!all.equal(chkPercT, wideIndivPercT))
  stop("Something different between the two sets of percentages.")
if (nrow(setdiff(chkPercT, wideIndivPercT)) != 0)
  stop("Extra observations were created when working with indivT")
if (nrow(setdiff(wideIndivPercT, chkPercT)) != 0)
  stop("More observations were are in the original worksheet than created when working with indivT")


rm(chkPercT)
## #######################
## ##################################################





## ##################################################
## Make other adjustments to the dataset so that it's easier to use.

## Remove the Bacteria row (last row), since it is just the totals of
## the taxa.  Remove the counts associated with unclassifed taxa.
## Also, include accum. degree days in the tibble.
indivT <- indivT %>%
  filter(!(origName %in% c("Bacteria", "Unclassified"))) %>%
  left_join(timeDF, by="days")


## Make a new, more readable taxa column.
## Column names with open brackets (e.g. "[Tissierellaceae]") causes
## problems for functions expecting traditional data frame column
## names.
indivT$taxa <- gsub(indivT$origName, pattern="\\[", replacement="")
indivT$taxa <- gsub(indivT$taxa, pattern="]", replacement="")
## Column names with dashes can likewise be a problem, so I replace
## dashes with underscores.
indivT$taxa <- gsub(indivT$taxa, pattern="-", replacement="_")
## Remova the taxonName column from the tibble to avoid confusion with
## the next taxa column.
indivT <- indivT %>% select(-origName)
## ##################################################




## ##################################################
## For use in graphs and in calculating percentages later, we need counts
## of total taxa types by:
##   Each pig and each day
##   Each day (all pigs combined)

## Total taxa counts by day (all pigs combined).
ctsByDayT <- indivT %>%
  group_by(days, degdays) %>%
  summarize(totals = sum(counts))
## ##################################################




## ##################################################
## Some taxa don't occur frequently.  We identify these, classify them
## as "rare", and then reformat the data accordingly.  Pechal et al
## (2013) call "rare" any taxa with <3% of relative abundance.  Also,
## we exclude all the "unclassified" counts.

## Find taxa which make up >= 3% of the total counts, over all days
## and subjects.  These are the "not rare" taxa, which I'll call
## "common".  Don't include the "unclassified" category of taxa, as
## unclassified taxa counts/fractions won't be used in the individual
## models.  (Note: Family-level taxa are unclassified much more often
## than 3%, whether you look at it by totals or by individual
## subj/day.)
commonByTotalV <- unlist(indivT %>%
    group_by(taxa) %>%
    summarize(taxatotal = sum(counts)) %>%
    mutate(frac = taxatotal/sum(taxatotal)) %>%
    filter(frac >= 0.03) %>%
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
## "unclassified" category of taxa.
commonByDaySubjV <- unlist(indivT %>%
                           left_join(ctsByDaySubjT) %>%
                           mutate(fracByDaySubj = counts/totals) %>%
                           group_by(taxa) %>%
                           summarize(maxFracByDaySubj = max(fracByDaySubj)) %>%
                           filter(maxFracByDaySubj >= 0.03) %>%
                           select(taxa)
                           )
## See a list of maximum taxa percentages sorted in descending order:
## data.frame(indivT %>%
##            left_join(ctsByDaySubjT) %>%
##            mutate(fracByDaySubj = counts/totals) %>%
##            group_by(taxa) %>%
##            summarize(maxFracByDaySubj = max(fracByDaySubj)) %>%
##            arrange(desc(maxFracByDaySubj))
##            )




## Yet another way to calculate the percentages represented by each
## taxa is to calculate the percentages by day (over all subjects).
## As above, we don't include the "unclassified" category of taxa.
commonByDayV <- unlist(indivT %>%
                       group_by(days, taxa) %>%
                       summarize(ctsByDayTaxa = sum(counts)) %>%
                       left_join(ctsByDayT) %>%
                       mutate(fracByDay = ctsByDayTaxa/totals) %>%
                       group_by(taxa) %>%
                       summarize(maxFracByDay = max(fracByDay)) %>%
                       filter(maxFracByDay >= 0.03) %>%
                       select(taxa)
                       )
## See a list of maximum taxa by day percentages sorted in descending
## order:
## data.frame(indivT %>%
##   group_by(days, taxa) %>%
##   summarize(ctsByDayTaxa = sum(counts)) %>%
##   left_join(ctsByDayT) %>%
##   mutate(fracByDay = ctsByDayTaxa/totals) %>%
##   group_by(taxa) %>%
##   summarize(maxFracByDay = max(fracByDay)) %>%
##   arrange(desc(maxFracByDay))
##   )
## ##################################################




## ##################################################
## Figure out which taxa should be treated as "rare".

## NOTE: Lists of common family-level taxa ARE different depending on
## how we calculate them!  It looks like the Pechal et al paper used the fractions calculated by day, so we go with that.
commonTaxaNamesV <- commonByDayV

rm(commonByDaySubjV, commonByDayV, commonByTotalV)



## Rename taxa that occur less than 3% of the time to "rare".  Then,
## sum the counts for these into one row.
commontaxaT <- indivT
commontaxaT[!(commontaxaT$taxa %in% commonTaxaNamesV), "taxa"] <- "Rare"
commontaxaT <- commontaxaT %>%
  group_by(days, degdays, subj, taxa) %>%
  summarize(counts = sum(counts))
rm(commonTaxaNamesV)
## WORKING HERE - Need to calculate the percentages by day and subject. ###
## ##################################################





## ##################################################
## Use the above total counts to find the fraction of the total
## represented by each taxa, whether for each pig/day or for each day
## (pigs combined).
indivT <- indivT %>%
  left_join(ctsByDaySubjT) %>%
  mutate(fracByDaySubj=counts/totals) %>%
  select(-totals)


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
