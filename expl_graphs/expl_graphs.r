library("tidyverse")
library("colorspace")


## ##################################################
## Are we dealing with phlya, orders, or families?
taxalevel <- "families"

## Read in cleaned-up phyla, orders, or families taxa.
taxaT <- read_csv(paste0("../initial_code/", taxalevel, "_massaged.csv"))
## ##################################################



## ##################################################
## Graphics illustrating the variability in counts.

## Using raw counts, we can see the large variability among individuals.
## For each bacteria, plot counts for each pig vs. day number.
ggplot(taxaT, aes(days, counts)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa) +
  xlab("Days") +
  ylab("Counts by degree day and subject") +
  labs(color="Subject")
## Save to a PDF file.
ggsave(paste0(taxalevel, "_scatter_counts_by_day_bacteria.pdf"), width=6, height=3.5, units="in")


## Make the stacked bar chart using raw counts by (all) degree days.
## ggplot(taxaT, aes(degdays)) +
##   geom_bar(aes(weight=counts, fill=taxa)) +
##   xlab("Degree days") +
##   ylab("Counts by degree day and taxa, combined over subjects")
## ##################################################



## ##################################################
## Look just at the first few days.  

## Subset to just the first few days.
days5T <- subset(taxaT, days <= 5)
## Add "Day" prefix to make graphs easier to read.
days5T$days <- paste0("Day ", days5T$days)

## Assess variability among pigs on the first few days, in terms of
## counts.
ggplot(days5T, aes(x=subj, y=counts, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days) +
  theme(axis.text.x = element_text(angle=90, hjust=0)) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(x="Subjects", y="Counts")

## Assess variability among pigs on the first few days, in terms of
## fractions.
ggplot(days5T, aes(x=subj, y=fracBySubjDay, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days) +
  theme(axis.text.x = element_text(angle=90, hjust=0)) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(x="Subjects", y="Composition fraction")
ggsave(paste0(taxalevel, "_first5_frac_bars_by_day_indiv.pdf"), width=6, height=3.5, units="in")


## Assess taxa variability across the first few days, averaged across
## pigs, in terms of fractions.  Compare this with Pechal et al,
## Fig. 1a.
avgdays5T <- taxaT %>%
  filter(days <= 5) %>%
  group_by(days, taxa) %>%
  summarize(avgFracByDay=mean(fracBySubjDay))
ggplot(avgdays5T, aes(x=days, y=avgFracByDay, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  labs(x="Days", y="Composition fraction")
ggsave(paste0(taxalevel, "_first5_avgfrac_bars_by_day.pdf"), width=3.5, height=3, units="in")

rm(days5T, avgdays5T)
## ##########################
## ##################################################




## ##################################################
## Scatterplot of percentage composition by degree day, with each
## bacteria a different color.  No distinguishing between individuals.

ggplot(taxaT %>% filter(taxa!="Rare"),
      aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=taxa)) +
  labs(x="Degree days", y="Composition fraction") +
  labs(color="Subject")
ggsave(paste0(taxalevel, "_scatter_frac_by_degday.pdf"), width=7, height=4, units="in")
## ##################################################



## ##################################################
## For each bacteria, scatterplot of fractions vs. accum. degree days
## (for each pig).

## For phylum taxa:
ggplot(taxaT %>% filter(taxa!="Rare"),
      aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa) +
  labs(x="Degree days", y="Composition fraction") +
  labs(color="Subject")
## ggsave(paste0(taxalevel, "_scatter_frac_by_degday_by_taxa.pdf"), width=6, height=4, units="in")
## ##################################################





## ########## STOPPED EDITING HERE! ##########




## ##################################################
## Bar charts using fractions associated with each taxa, rather than
## raw counts.

## Summarize the counts by day and by day-taxa, calculate percentages
## of each taxa per day (across all pigs).
bydayT <- taxaT %>%
  group_by(days, degdays, taxa) %>%
  summarize(ctsByTaxaDay = sum(counts)) %>%
  left_join(taxaT %>%
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
ggplot(subset(taxaT, days <= 5)),
       aes(x=subj, y=fracByDaySubj, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days)

ggplot(subset(taxaT, days <= 5),
       aes(x=subj, y=fracByDaySubj, fill=taxa)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(~days)
## ##################################################
