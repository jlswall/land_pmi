library("tidyverse")

## ##################################################
## Read the eukaryote and bacteria data which we previously organized
## separately.  Add a column to each to represent the type ("euk" or
## "bac"), so that when we combine, we can remember which taxa were of
## which type.
eukT <- read_csv("../../eukaryote_data/orders_massaged.csv")
eukT$type <- "euk"
bacT <- read_csv("../../revise_algorithm/orders_massaged.csv")
bacT$type <- "bac"


## Combine these two tibbles, renaming the "taxon" column in the
## eukaryote data to "taxa" (to match the way it's named in the
## bacteria tibble).
comboT <- rbind(eukT %>% rename(taxa=taxon), bacT)


## For the bacteria data, we did not have data for
## - subject A1 on days 7 and 9
## - subject A4 on day 7
## We did have data for these subjects and days in the eukaryote data,
## but because of the missing bacteria data, we remove these rows from
## the tibble.
comboT <- comboT %>% filter(! ((subj=="A1") & (days==7 | days==9)) )
comboT <- comboT %>% filter( !((subj=="A4") & (days==7)) )

## For the bacteria data, we had problems with the data from
## - subject A3 on day 40
## We did have data for these subjects and days in the eukaryote data,
## but because of the missing bacteria data, we remove these rows from
## the tibble.
comboT <- comboT %>% filter( !((subj=="A3") & (days==40)) )
## ##################################################



## ##################################################
## Add percentages by subj/day to the commontaxaT table.

## Check that the fractions add up to 1 on each day and subject for
## each type.
unique(
    unlist(comboT %>%
           group_by(type, days, subj) %>%
           summarize(sumFracBySubjDay = sum(fracBySubjDay)) %>%
           ungroup() %>%
           select(sumFracBySubjDay))
    )
## ##################################################



## ##################################################
## Save the tibble to a file for use in separate code
## for graphing and analysis.

write.csv(comboT, file="orders_massaged.csv", row.names=FALSE)
## ##################################################
