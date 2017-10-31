library("tidyverse")
library("colorspace")


## Read in cleaned-up phyla taxa.
phylumT <- read_csv("phyla_massaged.csv")
## Read in cleaned-up families taxa.
familyT <- read_csv("families_massaged.csv")



## ##################################################
## Graphics illustrating the variability in counts.

## Using raw counts, we can see the large variability among individuals.
## For each bacteria, plot counts for each pig vs. day number.
ggplot(phylumT, aes(days, counts)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa) +
  xlab("Days") +
  ylab("Counts by degree day and subject") +
  labs(color="Subject")
## Save to a PDF file.
ggsave("phyla_scatter_counts_by_day_bacteria.pdf", width=6, height=3.5, units="in")


## Make the stacked bar chart using raw counts by (all) degree days.
## ggplot(phylumT, aes(degdays)) +
##   geom_bar(aes(weight=counts, fill=taxa)) +
##   xlab("Degree days") +
##   ylab("Counts by degree day and taxa, combined over subjects")
## ##################################################



## ##################################################
## Look just at the first few days.  

## ##########################
## For phyla.

## counts and fractions.
phy5T <- subset(phylumT, days <= 5)
## Add "Day" prefix to make graphs easier to read.
phy5T$days <- paste0("Day ", phy5T$days)

## Assess variability among pigs on the first few days, in terms of
## counts.
## ggplot(phy5T, aes(x=subj, y=counts, fill=taxa)) +
##   geom_bar(stat="identity", position="stack") +
##   facet_grid(~days) +
##   theme(axis.text.x = element_text(angle=90, hjust=0)) +
##   scale_y_continuous(expand=c(0, 0)) +
##   labs(x="Subjects", y="Counts")

## Assess variability among pigs on the first few days, in terms of
## fractions.
ggplot(phy5T, aes(x=subj, y=fracBySubjDay, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days) +
  theme(axis.text.x = element_text(angle=90, hjust=0)) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(x="Subjects", y="Composition fraction")
ggsave("phyla_first5_frac_bars_by_day_indiv.pdf", width=6, height=3.5, units="in")


## Assess taxa variability across the first few days, averaged across
## pigs, in terms of fractions.  Compare this with Pechal et al,
## Fig. 1a.
avgPhy5T <- phylumT %>%
  filter(days <= 5) %>%
  group_by(days, taxa) %>%
  summarize(avgFracByDay=mean(fracBySubjDay))
ggplot(avgPhy5T, aes(x=days, y=avgFracByDay, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  labs(x="Days", y="Composition fraction")
ggsave("phyla_first5_avgfrac_bars_by_day.pdf", width=3.5, height=3, units="in")

rm(phy5T, avgPhy5T)
## ##########################


## ##########################
## For families.

## counts and fractions.
fam5T <- subset(familyT, days <= 5)
## Add "Day" prefix to make graphs easier to read.
fam5T$days <- paste0("Day ", fam5T$days)

## Assess variability among pigs on the first few days, in terms of
## counts.
## ggplot(fam5T, aes(x=subj, y=counts, fill=taxa)) +
##   geom_bar(stat="identity", position="stack") +
##   facet_grid(~days) +
##   theme(axis.text.x = element_text(angle=90, hjust=0)) +
##   scale_y_continuous(expand=c(0, 0)) +
##   labs(x="Subjects", y="Counts")

## Assess variability among pigs on the first few days, in terms of
## fractions.
## ggplot(fam5T, aes(x=subj, y=fracBySubjDay, fill=taxa)) +
##   geom_bar(stat="identity", position="stack") +
##   facet_grid(~days) +
##   theme(axis.text.x = element_text(angle=90, hjust=0)) +
##   scale_y_continuous(expand=c(0, 0)) +
##   labs(x="Subjects", y="Composition fraction")
## ggsave("families_first5_frac_bars_by_day_indiv.pdf", width=6, height=3.5, units="in")


## Assess taxa variability across the first few days, averaged across
## pigs, in terms of fractions.  Compare this with Pechal et al,
## Fig. 1b.
avgFam5T <- familyT %>%
  filter(days <= 5) %>%
  group_by(days, taxa) %>%
  summarize(avgFracByDay=mean(fracBySubjDay))
ggplot(avgFam5T, aes(x=days, y=avgFracByDay, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  labs(x="Days", y="Composition fraction")
ggsave("families_first5_avgfrac_bars_by_day.pdf", width=4.5, height=5.5, units="in")

rm(fam5T, avgFam5T)
## ##########################
## ##################################################



## ########## STOPPED EDITING HERE! ##########




## ##################################################
## Bar charts using fractions associated with each taxa, rather than
## raw counts.

## Summarize the counts by day and by day-taxa, calculate percentages
## of each taxa per day (across all pigs).
bydayT <- phylumT %>%
  group_by(days, degdays, taxa) %>%
  summarize(ctsByTaxaDay = sum(counts)) %>%
  left_join(phylumT %>%
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

## ####### WORKING HERE! 



## Bar charts with counts.
ggplot(subset(phylumT, days <= 5)),
       aes(x=subj, y=fracByDaySubj, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days)

ggplot(subset(phylumT, days <= 5),
       aes(x=subj, y=fracByDaySubj, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days)
## ##################################################
