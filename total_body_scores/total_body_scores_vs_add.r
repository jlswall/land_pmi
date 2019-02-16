library("tidyverse")
library("ggplot2")
library("readxl")
library("cowplot")

## ##################################################
## Read in data from Excel.

## First, read in ADDs.
fileNm <- "total_body_scores.xlsx"
listADDs <- as.vector(as.matrix(read_excel(path=fileNm, range="D1:S1", col_names=F)))

## Read in the rows with the TBS for the 6 cadavers.
rawAllT <- read_excel(path=fileNm, range="A22:S27", col_names=F)
## Drop column B which contains just the string "TBS".  Drop column C,
## which is all 0s (TBS not meaningful because this is the day on
## which cadavers were placed).
wideT <- rawAllT[,-c(2, 3)]
colnames(wideT) <- c("subj", paste("ADD", listADDs, sep="_"))
rm(fileNm, listADDs)

## Now, we need to go from wide to long format.
intermedT <- gather(wideT, ADD_degdays, tbs, -subj)
## In the degdays column, we need to remove the "ADD_" prefix, so that
## we can get the numeric results (accumulated degree days).
tbsADDT <- separate(intermedT, ADD_degdays, c("add", "degdays"), sep="_", convert=T)
## Remove the column that contains only the old "ADD" prefix.
tbsADDT <- tbsADDT %>% select(-add)
rm(intermedT)
## ##################################################



## ##################################################
## Make scatterplot of TBS vs. log10(ADD), with regression line,
## regression equation, and R squared superimposed.

## For help with setting up the expression for the regression equation, see
## https://stackoverflow.com/questions/7549694/adding-regression-line-equation-and-r2-on-graph
linLM <- lm(tbs ~ log10(degdays), data=tbsADDT)
l <- list(a = round(as.vector(coef(linLM))[1], 4),
          b = round(as.vector(coef(linLM))[2], 4))
eq <- as.character(as.expression(substitute(paste(hat(y) == a + b, "x"), l)))
Rsq <- round(summary(linLM)$r.squared, 4)
lPanel <- ggplot(tbsADDT, aes(x=log10(degdays), y=tbs)) +
  geom_point() +
  stat_smooth(method="lm", se=FALSE) +
  annotate("text", x=2.2, y=7.0, hjust=0, label=eq, parse=TRUE, size=3) +
  annotate("text", x=2.2, y=5.5, hjust=0,
           label=paste("R^2  ==", Rsq), parse=T, size=3) +
  theme_bw() +
  theme(plot.margin=unit(rep(0.15, 4), "in")) +
  labs(x="Log 10 Accumulated Degree Days", y="Total Body Score")

rm(linLM, l, eq, Rsq)
## ##################################################



## ##################################################
## Make scatterplot of TBS vs. ADD, with regression line,
## regression equation, and R squared superimposed.

## For help with setting up the expression for the regression equation, see
## https://stackoverflow.com/questions/7549694/adding-regression-line-equation-and-r2-on-graph
curveLM <- lm(tbs ~ degdays + I(degdays^2) + I(degdays^3), data=tbsADDT)
## l <- list(a = round(as.vector(coef(curveLM))[1], 4),
##           b1 = round(as.vector(coef(curveLM))[2], 4),
##           b2 = round(as.vector(coef(curveLM))[3], 4),
##           b3 = round(as.vector(coef(curveLM))[4], 4))
paramLst <- as.list(as.numeric(format(c(
    as.vector(coef(curveLM))[1],
    as.vector(coef(curveLM))[2],
    as.vector(coef(curveLM))[3],
    as.vector(coef(curveLM))[4]), digits=4)))
names(paramLst) <- c("a", "b1", "b2", "b3")
eq <- as.character(as.expression(substitute(paste(hat(y) == a + b1, "x", b2, x^2 + b3, x^3), paramLst)))
## eq <- as.character(as.expression(substitute(paste(hat(y) == a + b1, "x", b2, x^2 +, b3, x^3), paramLst)))
Rsq <- round(summary(curveLM)$r.squared, 4)
rPanel <- ggplot(tbsADDT, aes(x=degdays, y=tbs)) +
  geom_point() +
  stat_smooth(method="lm", formula="y ~ x + I(x^2) + I(x^3)", se=FALSE) +
  annotate("text", x=168, y=7, hjust=0, label=eq, parse=TRUE, size=3) +
  annotate("text", x=168, y=5.5, hjust=0,
           label=paste("R^2  ==", Rsq), parse=T, size=3) +
  theme_bw() +
  theme(plot.margin=unit(rep(0.15, 4), "in")) +
  labs(x="Accumulated Degree Days", y="Total Body Score")
## ##################################################


## ##################################################
plot_grid(lPanel, rPanel, labels=c("a", "b"))
ggsave(file="tbs_vs_add_2panels.pdf", height=4, width=8, units="in")
## ##################################################



## ##################################################
## Code I was hoping to use to combine our ADDs with the total body
## scores doesn't work, because there are 19 days of total body scores
## and only 18 ADDs.

## ## Get the days and accumulated degree days from the other data we've
## ## already read in.
## fileNm <- "families_massaged.csv"
## rawdaysT <- read_csv(file=fileNm)
## daysADDT <- unique(rawdaysT %>% select(days, degdays))
## rm(rawdaysT, fileNm)


## rawAllT <- read_excel(path=fileNm, range="A22:S27", col_names=F)
## ## Drop the column which contains just the string "TBS".
## wideT <- rawAllT[,-2]
## colnames(wideT) <- c("subj", paste("ADD.", daysADDT$degdays, sep=""))
## rm(fileNm)

## ## Column names are unique in this sheet.
## sum(duplicated(colnames(rawAllT)))

## ## Print out our mismatched ADDs and these, side by side.
## write.table(file="mismatched_ADD_in_tbs_file.txt", cbind(c(daysADDT$degdays, NA), c(NA, listADDs)), sep="\t", row.names=F, col.names=F)
