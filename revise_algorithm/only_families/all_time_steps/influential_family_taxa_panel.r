library("tidyverse")

## These are the top 5 families (in terms of IncNodePurity) from the
## model with the stricter inclusion criterion (hit_1perc_twice):
topChoices <- c("Peptostreptococcaceae", "Microbacteriaceae", "Turicibacteraceae", "Oxalobacteraceae", "Bifidobacteriaceae")
## We also add Listeriaceae, which was important in our original
## analysis, but has been excluded under the strict criterion.
topChoices <- c(topChoices, "Listeriaceae")


## Using the dataset with the weird observations from subj. A3, day 40.
indivT <- read_csv("../../families_massaged.csv")

chooseT <- indivT %>%
  filter(taxa %in% topChoices)
chooseT$taxa <- factor(chooseT$taxa, levels=topChoices)

ggplot(chooseT, aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  labs(x="Degree days", y="Fraction") +
  facet_wrap(~taxa, scales="free_y")  ## Allow diff. y-scales across panels.
  ## facet_wrap(~taxa)  ## Keep y-scales same across panels.
ggsave("influential_family_taxa_panel.pdf", width=7, height=4, units="in")
