library("tidyverse")

## These are the top 6 orders (in terms of IncNodePurity) from the
## model with the stricter inclusion criterion (hit_1perc_twice):
topChoices <- c("Turicibacterales", "Lactobacillales", "Burkholderiales", "JG30_KF_CM45", "Sphingobacteriales", "Rhizobiales")


## Using the dataset without the weird observations from subj. A3, day 40.
indivT <- read_csv("../../orders_massaged.csv")

chooseT <- indivT %>%
  filter(taxa %in% topChoices)
chooseT$taxa <- factor(chooseT$taxa, levels=topChoices)

ggplot(chooseT, aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  labs(x="Degree days", y="Fraction") +
  facet_wrap(~taxa, scales="free_y")  ## Allow diff. y-scales across panels.
  ## facet_wrap(~taxa)  ## Keep y-scales same across panels.
ggsave("influential_order_taxa_panel.pdf", width=7, height=4, units="in")
