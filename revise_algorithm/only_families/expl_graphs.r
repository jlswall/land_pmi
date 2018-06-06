library("tidyverse")

## Using the dataset with the weird observations from subj. A3, day 40.
indivT <- read_csv("../with_weird_subjA3day40_families_massaged.csv")

freqtaxaT <- indivT %>%
  group_by(taxa) %>%
  summarize(maxFracBySubjDay=max(fracBySubjDay)) %>%
  arrange(desc(maxFracBySubjDay))


graph1T <- freqtaxaT %>%
  filter(maxFracBySubjDay >= 0.10) %>%
  left_join(indivT) %>%
  select(-maxFracBySubjDay)
ggplot(graph1T, aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa)
ggsave("family_taxa_w_maxfrac_at_least_10perc.pdf")


graph2T <- freqtaxaT %>%
  filter((maxFracBySubjDay < 0.10) & (maxFracBySubjDay >= 0.05)) %>%
  left_join(indivT) %>%
  select(-maxFracBySubjDay)
ggplot(graph2T, aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa)
ggsave("family_taxa_w_maxfrac_5-10perc.pdf")


graph3T <- freqtaxaT %>%
  filter(maxFracBySubjDay < 0.05) %>%
  left_join(indivT) %>%
  select(-maxFracBySubjDay)
ggplot(graph3T, aes(degdays, fracBySubjDay)) +
  geom_point(aes(color=subj)) +
  facet_wrap(~taxa)
ggsave("family_taxa_w_maxfrac_less_5perc.pdf")
