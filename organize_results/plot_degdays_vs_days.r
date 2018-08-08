library("figdim")

## Read in data from CSV file.
daysDF <- read.csv("degdays_vs_days.csv")

## Make a column indicating whether the data point falls during the
## first 15 days.
daysDF$isFirst15Days <- ifelse(daysDF$Day<=15, 1, 0)


init.fig.dimen(file="degdays_vs_days.pdf", f.type="pdf", width=3.75, height=4.25, mai=c(0.6, 0.6, 0.05, 0.05))
plot(daysDF$Day, daysDF$Degree.day, type="n", xlab="Days post mortem", ylab="Accumulated degree days post mortem")
## The estimated regression equation is:
## Estimated Degree.day = 10.7 + (27.9 * Day)
## lines(x=c(17, 17), y=c(-200, 485.2), col="gray")
## lines(x=c(-20, 17), y=c(485.2, 485.2), col="gray")
## Draw in points for the first 15 days.
with(subset(daysDF, isFirst15Days==1), points(Day, Degree.day))
## Draw in points for the remaining days.
with(subset(daysDF, isFirst15Days==0), points(Day, Degree.day, pch=19))
## Give legend to distinguish the points.
legend("bottomright", legend=c("within first 15 days", "after first 15 days"), pch=c(1, 19), cex=0.8, pt.cex=0.8)
dev.off()

## Correlation between days and degree days for the whole time period.
with(daysDF, cor(Day, Degree.day))
## 0.9998702
## Correlation between days and degree days for just the first 15 days.
with(subset(daysDF, isFirst15Days==1), cor(Day, Degree.day))
## 0.9999598
