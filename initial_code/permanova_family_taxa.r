library("vegan")
library("colorspace")



## ##################################################
## Treat each pig as a "site"/"community" and each day as a time step.
## The repsonse is multivariate because there are counts for many
## different phyla.

## We need to get the data into "community matrix format", which is
## basically a wide format.
wideT <- commontaxaT %>%
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
