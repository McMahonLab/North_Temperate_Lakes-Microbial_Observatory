#LSA metric code

#1. Read in LSA output
#2. Measure the connectitivity of each OTU
#   2a. do this in binary - number of edges
#   2b. do this more quantitatively - strength and number of edges
#3. Compute a metric of connectivity for each sample over time
#   use sum(abun x number of edges/strength of edges)

setwd("C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/Network_analysis")
TBH.network <- read.table(file = "TBH_network_25Jan16.txt", header=T)

#For each otu, calculate the number of edges.
#This is equivalent to the number of times it appears in the file

TBH.edges <- table(c(as.character(TBH.network$index1), as.character(TBH.network$index2)))

#Quantify the strength of these connections - take average correlation * number of edges

TBH.corr <- c()
for(i in 1:length(TBH.edges)){
  hits <- which(TBH.network$index1 == names(TBH.edges)[i] | TBH.network$index2 == names(TBH.edges)[i])
  TBH.corr[i] <- mean(TBH.network$LSA[hits])
}

TBH.quant <- TBH.edges * TBH.corr

#For each TBH sample, take the abundance of each OTU * the connectivity metric
library(OTUtable)
data(otu_table)
TBH <- bog_subset("TBH", otu_table)

nodes <- match(names(TBH.edges), rownames(TBH))
TBH.conn <- TBH[nodes,]
metric1 <- colSums(sweep(TBH.conn, 1, TBH.edges, "*"))
metric2 <- colSums(sweep(TBH.conn, 1, TBH.quant, "*"))

TBH.dates <- extract_date(colnames(TBH))

plot(TBH.dates, metric2)

#Both metrics show the same trends. Stick with metric1.

plot.conn <- data.frame(TBH.dates, metric1)
library(ggplot2)
ggplot(data=plot.conn, aes(x=TBH.dates, y=metric1)) + geom_line()
