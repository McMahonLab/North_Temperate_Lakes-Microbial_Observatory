# Script to generate figures for NTL-MO manuscript
# Plots are saved as pdfs. For many figures, further processing is performed in Adobe Illustrator
###########
# Input your path to the North_Temperate_Lakes-Microbial_Observatory/Figures folder in the GitHub repo, or whereever you would like the figures to be saved
path2repo <- "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/"

# Load packages
library(OTUtable)     # Contains data and functions for analysis of NTL-MO OTU table
library(exactRankTests)# Calculates Wilcoxon Rank Significance for ties on richness betwen lakes
library(vegan)        # Used for Bray-Curtis
library(ggplot2)      # Used for plotting
library(reshape2)     # Used to format metadata
library(grid)         # Used in multiplot()
library(ggdendro)     # Used to plot hierarchical clustering trees
library(indicspecies) # Performs indicator taxa analysis

# Load data from OTUtable
data(otu_table)
data(taxonomy)
data(metadata)

# Generate tables at the clade and phylum level
clade_table <- combine_otus("Clade", otu_table, taxonomy)
phylum_table <- combine_otus("Phylum", otu_table, taxonomy)

# Reduce names to last known taxonmic information
clade_table <- reduce_names(clade_table)
phylum_table <- reduce_names(phylum_table)
clade_table07 <- year_subset("07", clade_table)

# Set the multiplot function (from user on Stack Overflow). Used to output multiple plots in one pdf
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])   
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

###############

##############
# Figure 1
# List lake categories
lakes <- c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA")

# Split epilimnion and hypolimnion into separate tables
epilimnia <- bog_subset("..E", otu_table)
hypolimnia <- bog_subset("..H", otu_table)

# Calculate observed richness
epi.obs <- apply(epilimnia, 2, obs_richness)
hypo.obs <- apply(hypolimnia, 2, obs_richness)

# Extract sampling location from sample names
epi.lakes <- substr(names(epi.obs), start = 1, stop = 2)
hypo.lakes <- substr(names(hypo.obs), start = 1, stop = 2)

# Make dataframe for plotting
epi.data <- data.frame(epi.lakes, epi.obs)
hypo.data <- data.frame(hypo.lakes, hypo.obs)
# Order factors
epi.data$epi.lakes <- ordered(epi.data$epi.lakes, levels = lakes)
hypo.data$hypo.lakes <- ordered(hypo.data$hypo.lakes, levels = lakes)

# plot 1A
pdf(file = paste(path2repo, "epi_boxplot.pdf", sep = ""), width = 3.3125, height = 3)
ggplot(data = epi.data, aes(y = epi.obs, x = epi.lakes, fill = epi.lakes)) + geom_boxplot() + labs(y = "Observed Richness", x = NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), legend.position = "none") + scale_fill_brewer(palette = "Set3")
dev.off()

# plot 1B
pdf(file = paste(path2repo, "hypo_boxplot.pdf", sep = ""), width = 3.3125, height = 3)
ggplot(data = hypo.data, aes(y = hypo.obs, x = hypo.lakes, fill = hypo.lakes)) + geom_boxplot() + labs(y="Observed Richness", x = NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title.y = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), legend.position="none") + scale_fill_brewer(palette = "Set3")
dev.off()

# Check significance using the Wilcoxon Rank Test on medians (not output as pdf, indicated as symbols in Illustrator)

lakeids <- levels(epi.data$epi.lakes)
Lake1 <- character(0)
Lake2 <- character(0)
Test <- character(0)
pvalue <- numeric(0)
Interpretation <- character(0)

for (id in lakeids) { 
  x = subset(epi.data, epi.data$epi.lakes == id)
  for (lake in lakeids) {
    y = subset(epi.data, epi.data$epi.lakes == lake)
    result <- wilcox.exact(x$epi.obs,y$epi.obs, alternative= "two.sided", conf.level = 0.95)
    Lake1 <- c(Lake1, id)
    Lake2 <- c(Lake2, lake)
    Test <- c(Test, result$method)
    pvalue <- c(pvalue, result$p.value)
    if (result$p.value < 0.05) {
      Interpretation <- c(Interpretation, "Reject Null (difference in median detected)")
    } else {
      Interpretation <- c(Interpretation, "Accept Null (difference in median undetected)")
    }  
  }
}

Wilcoxon.epi <- data.frame(Lake1, Lake2, Test, pvalue, Interpretation)
print(Wilcoxon.epi)

lakeids <- levels(hypo.data$hypo.lakes)
Lake1 <- character(0)
Lake2 <- character(0)
Test <- character(0)
pvalue <- numeric(0)
Interpretation <- character(0)

for (id in lakeids) { 
  x = subset(hypo.data, hypo.data$hypo.lakes == id)
  for (lake in lakeids) {
    y = subset(hypo.data, hypo.data$hypo.lakes == lake)
    result <- wilcox.exact(x$hypo.obs,y$hypo.obs, alternative= "two.sided", conf.level = 0.95)
    Lake1 <- c(Lake1, id)
    Lake2 <- c(Lake2, lake)
    Test <- c(Test, result$method)
    pvalue <- c(pvalue, result$p.value)
    if (result$p.value < 0.05) {
      Interpretation <- c(Interpretation, "Reject Null (difference in median detected)")
    } else {
      Interpretation <- c(Interpretation, "Accept Null (difference in median undetected)")
    }  
  }
}

Wilcoxon.hypo <- data.frame(Lake1, Lake2, Test, pvalue, Interpretation)
print(Wilcoxon.hypo)
###################
# Figure 2
# Trout Bog, 2007

# Identify mixing dates (less than 1 degree of temperature difference between 0.5 meters and maximum sampling depth)
metalakes <- substr(metadata$Sample_Name, start = 1, stop = 3)
metayears <- substr(metadata$Sample_Name, start = 9, stop = 10)
metaTBH <- metadata[which(metalakes == "TBH" & metayears == "07"), c(1,2,4)]
metaTBH <- dcast(metaTBH, Sample_Name ~ Depth, fun.aggregate = mean)
TBHmixes <- extract_date(metaTBH$Sample_Name[which(metaTBH$"0.5" - metaTBH$"7" < 1)])

# Make dataset of Trout Bog hypolimion samples from 2007
hypo <- bog_subset(paste("TBH", sep = ""), otu_table)
hypo <- year_subset("07", hypo)
# Calculate observed richness
hypo.rich <- apply(hypo, 2, obs_richness)
# Extract sampling date from sample names
hypo.date <- extract_date(colnames(hypo))
# Remove January samples - large gap distracts in plot, and winter samples are not considered in this study
hypo.rich <- hypo.rich[c(1:32, 35:80)]
hypo.date <- hypo.date[c(1:32, 35:80)]

# Make dataframe for plotting
TB_richness <- data.frame(hypo.date, hypo.rich)
colnames(TB_richness) <- c("date", "richness")

pdf(file = paste(path2repo, "richness_over_time1.pdf", sep = ""), width = 3.3125, height = 2.3125)
ggplot() + geom_line(data = TB_richness, aes(x = date, y = richness), size = 1) + labs(title = "Trout Bog", x = NULL, y = "Observed Richness") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks=element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title = element_text(size = 10, vjust = 2), axis.text.y = element_text(colour = "black", size=10), plot.title = element_text(size = 12, vjust = 2), legend.position = "none") + geom_point(data = TB_richness[match(TBHmixes, TB_richness$date), ], aes(x = date, y = richness), size = 2, colour = "red")
dev.off()


# Repeat with North Sparkling, 2008
metaNSH <- metadata[which(metalakes == "NSH" & metayears == "08"), c(1,2,4)]
metaNSH <- dcast(metaNSH, Sample_Name ~ Depth, fun.aggregate = mean)
NSHmixes <- extract_date(metaNSH$Sample_Name[which(metaNSH$"0.5" - metaNSH$"4" < 1)])

hypo <- bog_subset(paste("NSH", sep = ""), otu_table)
hypo <- year_subset("08", hypo)
hypo.rich <- apply(hypo, 2, obs_richness)
hypo.date <- extract_date(colnames(hypo))

NS_richness <- data.frame(hypo.date, hypo.rich)
colnames(NS_richness) <- c("date", "richness")

pdf(file = paste(path2repo, "richness_over_time2.pdf", sep = ""), width = 3.3125, height = 2.3125)
ggplot() + geom_line(data = NS_richness, aes(x = date, y = richness), size = 1.2) + labs(title = "North Sparkling Bog", x = NULL, y = "Observed Richness") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title = element_text(size = 10, vjust=2), axis.text.y = element_text(colour = "black", size = 10), plot.title = element_text(size=12, vjust = 2), legend.position = "none") + geom_point(data = NS_richness[match(NSHmixes, NS_richness$date), ], aes(x = date, y = richness), size = 2, colour = "red")
dev.off()


#############
# Figure 3 - Network analysis
all.network <- read.table(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Network_analysis/allsamples_network_28Jan16.txt", header = T)
all.edges <- table(c(as.character(all.network$index1), as.character(all.network$index2)))

# Get OTU abundance data
TBH <- bog_subset("TBH", otu_table)
NSH <- bog_subset("NSH", otu_table)
MAH <- bog_subset("MAH", otu_table)

# Calcuate connectivity metric for each sample by lake
# metric (for each sample) = sum(OTU abundance * # of connections to that OTU * average strength of connections to that OTU)
TBH.nodes <- match(names(all.edges), rownames(TBH))
TBH.conn <- TBH[TBH.nodes,]
TBH.edges <- all.edges[match(rownames(TBH.conn), names(all.edges))]

TBH.corr <- c()
for(i in 1:length(TBH.edges)){
  hits <- which(all.network$index1 == names(TBH.edges)[i] | all.network$index2 == names(TBH.edges)[i])
  TBH.corr[i] <- mean(all.network$LSA[hits])
}
TBH.quant <- TBH.edges * TBH.corr
TBH.metric <- colSums(sweep(TBH.conn, 1, TBH.quant, "*"))
TBH.dates <- extract_date(colnames(TBH))

NSH.nodes <- match(names(all.edges), rownames(NSH))
NSH.conn <- NSH[NSH.nodes,]
NSH.edges <- all.edges[match(rownames(NSH.conn), names(all.edges))]

NSH.corr <- c()
for(i in 1:length(NSH.edges)){
  hits <- which(all.network$index1 == names(NSH.edges)[i] | all.network$index2 == names(NSH.edges)[i])
  NSH.corr[i] <- mean(all.network$LSA[hits])
}
NSH.quant <- NSH.edges * NSH.corr
NSH.metric <- colSums(sweep(NSH.conn, 1, NSH.quant, "*"))
NSH.dates <- extract_date(colnames(NSH))

MAH.nodes <- match(names(all.edges), rownames(MAH))
MAH.conn <- MAH[MAH.nodes,]
MAH.edges <- all.edges[match(rownames(MAH.conn), names(all.edges))]

MAH.corr <- c()
for(i in 1:length(MAH.edges)){
  hits <- which(all.network$index1 == names(MAH.edges)[i] | all.network$index2 == names(MAH.edges)[i])
  MAH.corr[i] <- mean(all.network$LSA[hits])
}
MAH.quant <- MAH.edges * MAH.corr
MAH.metric <- colSums(sweep(MAH.conn, 1, MAH.quant, "*"))
MAH.dates <- extract_date(colnames(MAH))


all.metric <- c(TBH.metric, NSH.metric, MAH.metric)
all.dates <- c(TBH.dates, NSH.dates, MAH.dates)
lakekey <- c(rep("TBH", length(TBH.metric)), rep("NSH", length(NSH.metric)), rep("MAH", length(MAH.metric)))
plot.conn <- data.frame(lakekey, all.dates, all.metric)
colnames(plot.conn) <- c("Lake", "Date", "Connectivity")
pdf(file = paste(path2repo, "connectivity2007.pdf", sep = ""), width = 3.3125 * 2, height = 2.3)
ggplot(data = plot.conn, aes(x = Date, y = Connectivity, colour = Lake)) + geom_line(size = 1) + scale_y_log10() + coord_cartesian(xlim = extract_date(c("TBH15May07", "TBH18Nov07"))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title.x = element_text(size = 12, vjust = 0.3), axis.title.y = element_text(size=12, vjust=1.3), axis.text.y = element_text(colour = "black", size = 10))
dev.off()

pdf(file = paste(path2repo, "connectivity2008.pdf", sep = ""), width = 3.3125 * 2, height = 2.3)
ggplot(data = plot.conn, aes(x = Date, y = Connectivity, colour = Lake)) + geom_line(size = 1) + scale_y_log10() + coord_cartesian(xlim = extract_date(c("TBH15May08", "TBH18Nov08"))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title.x = element_text(size = 12, vjust=0.3), axis.title.y = element_text(size = 12, vjust = 1.3), axis.text.y = element_text(colour = "black", size = 10))
dev.off()


# Use linear model to show significance

# Make a variable of continous mixing
metalakes <- substr(metadata$Sample_Name, start=1, stop=3)
metayears <- substr(metadata$Sample_Name, start=9, stop=10)

metaTBH <- metadata[which(metalakes == "TBH"), c(1,2,3)]
metaTBH <- dcast(metaTBH, Sample_Name~Depth, fun.aggregate=mean)

metaNSH <- metadata[which(metalakes == "NSH"), c(1,2,3)]
metaNSH <- dcast(metaNSH, Sample_Name~Depth, fun.aggregate=mean)

metaMAH <- metadata[which(metalakes == "MAH"), c(1,2,3)]
metaMAH <- dcast(metaMAH, Sample_Name~Depth, fun.aggregate=mean)

TBHdates <- extract_date(metaTBH$Sample_Name)
NSHdates <- extract_date(metaNSH$Sample_Name)
MAHdates <- extract_date(metaMAH$Sample_Name)

cont.mixes <- c()
for(i in 1:dim(plot.conn)[1]){
  if(plot.conn$Lake[i] == "TBH"){
    sample <- metaTBH[which(TBHdates == plot.conn$Date[i]),]
    sample <- sample[2:length(sample)]
    sample <- as.numeric(sample[which(is.na(sample) == F)])
    cont.mixes[i] <- sample[3] - min(sample)
  }else if(plot.conn$Lake[i] == "NSH"){
    sample <- metaNSH[which(NSHdates == plot.conn$Date[i]),]
    sample <- sample[2:length(sample)]
    sample <- as.numeric(sample[which(is.na(sample) == F)])
    cont.mixes[i] <- sample[3] - min(sample)
  }else if(plot.conn$Lake[i] == "MAH"){
    sample <- metaMAH[which(MAHdates == plot.conn$Date[i]),]
    sample <- sample[2:length(sample)]
    sample <- as.numeric(sample[which(is.na(sample) == F)])
    cont.mixes[i] <- sample[3] - min(sample)
  }
}

plot.conn$Mixing.cont <- cont.mixes

# Add year and Julian Date variables
year <- substr(plot.conn$Date, start=1, stop=4)
plot.conn$Year<- factor(year, levels=c("2005", "2007", "2008", "2009"))

julian <- c()
for(i in 1:dim(plot.conn)[1]){
  if(year[i] == "2005"){
    julian[i] <- as.numeric(plot.conn$Date[i] - extract_date(c("TBH01Jun05")))
  }else if(year[i] == "2007"){
    julian[i] <- as.numeric(plot.conn$Date[i] - extract_date(c("TBH01Jun07")))
  }else if(year[i] == "2008"){
    julian[i] <- as.numeric(plot.conn$Date[i] - extract_date(c("TBH01Jun08")))
  }else if(year[i] == "2009"){
    julian[i] <- as.numeric(plot.conn$Date[i] - extract_date(c("TBH01Jun09")))
  }
}

plot.conn$DayNum <- julian

# Run linear model
# Setup: Taking log transformed connectivity (+1 to avoid issues with 0) is determined by the year, the lake, and the interaction of the Julian date and the continuous mixing variable ()
model <- lm(log(Connectivity+1) ~Year + Lake + Mixing.cont*DayNum, data=plot.conn[which(plot.conn$DayNum >= 0),])
# Run model without interaction to see main effects
#summary(lm(log(Connectivity+1) ~Year + Lake + Mixing.cont+DayNum, data=plot.conn[which(plot.conn$DayNum >= 0),]))
summary(model)

# Check the residuals
plot(resid(model) ~ fitted(model))

# The model tells me that day number effects connectivity differently depending on how strongly stratified the lake is
# Check the predicted values for how day number affects different lakes
plot(fitted(model)[which(plot.conn$DayNum >= 0 & plot.conn$Lake == "MAH")] ~ plot.conn$DayNum[which(plot.conn$DayNum >= 0 & plot.conn$Lake == "MAH")])
# Little to no trend of day num in Mary
plot(plot.conn$DayNum[which(plot.conn$Lake == "TBH" & plot.conn$Mixing.cont > 3 )], plot.conn$Connectivity[which(plot.conn$Lake == "TBH" & plot.conn$Mixing.cont > 3 )])
# Slight increase over time when lake is stratified
plot(plot.conn$DayNum[which(plot.conn$Lake == "TBH" & plot.conn$Mixing.cont < 3 )], plot.conn$Connectivity[which(plot.conn$Lake == "TBH" & plot.conn$Mixing.cont < 3 )])
# Seems to be either high or low. Likely some cool epilimnion samples right before stratification are contributing high connectivity, but the lake would still be stratified (but just barely)

plot(plot.conn$DayNum[which(plot.conn$Lake == "NSH" & plot.conn$Mixing.cont > 3 )], plot.conn$Connectivity[which(plot.conn$Lake == "NSH" & plot.conn$Mixing.cont > 3 )])
#Little to no trend over time

plot(plot.conn$DayNum[which(plot.conn$Lake == "NSH" & plot.conn$Mixing.cont < 3 )], plot.conn$Connectivity[which(plot.conn$Lake == "NSH" & plot.conn$Mixing.cont < 3 )])
# Mixed samples almost uniformly low

# Conclusions:
# a) Mary has higher connectivity than the other lakes
# b) Connectivity is different if the lake is mixed vs stratified (observationally, lower when mixed)
# c) The effect of date on connectivity is different if the lake is mixed vs stratified (observationally, increase over time when stratified, no trend when mixed)


#############
# Figure 4 - dendograms

TBH <- bog_subset("TBH", otu_table)
TBH <- TBH[which(rowSums(TBH) > 100), ]
input <- remove_reps(TBH)
month <- substr(colnames(input), start = 6, stop = 8)
input <- input[, which(month == "MAY" | month == "JUN" | month == "AUG" | month == "JUL" | month == "SEP")]
colnames(input) <- substr(colnames(input), start = 4, stop = 10)
h <- hclust(vegdist(t(input), method = "bray"), method = "ward.D2")

ddata <- dendro_data(h, type = "rectangle")
ddata$year <- factor(substr(ddata$label$label, start = 6, stop = 7), levels = c("05", "07", "08", "09"))

pdf(file = paste(path2repo, "TBH_dendrogram_abundancecutoff.pdf", sep = ""), width = 3.3125 * 2, height = 2.3)
ggplot(segment(ddata)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", angle = 90)) + geom_text(data = label(ddata), aes(x = x, y = y, label = label, colour = ddata$year), vjust = 0.5, hjust = 1.1, size = 2, angle = 90)   + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"), legend.position = "none", axis.text.y = element_text(angle = 0)) + scale_color_manual(values = c("darkorange3", "chartreuse4", "darkslategrey", "darkred"))  + coord_cartesian(ylim = c(-1.5, 2))
dev.off()

MAH <- bog_subset("MAH", otu_table)
MAH <- MAH[which(rowSums(MAH) > 100), ]
input <- remove_reps(MAH)
month <- substr(colnames(input), start = 6, stop = 8)
input <- input[,which(month == "MAY" | month == "JUN" | month == "AUG" | month == "JUL" | month == "SEP")]
colnames(input) <- substr(colnames(input), start = 4, stop = 10)
h <- hclust(vegdist(t(input), method = "bray"), method = "ward.D2")

ddata <- dendro_data(h, type = "rectangle")
ddata$year <- factor(substr(ddata$label$label, start = 6, stop = 7), levels = c("05", "07", "08", "09"))

pdf(file = paste(path2repo, "MAH_dendrogram_abundancecutoff.pdf", sep = ""), width = 3.3125 * 2, height = 2.3)
ggplot(segment(ddata)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", angle = 90)) + geom_text(data = label(ddata), aes(x = x, y = y, label = label, colour = ddata$year), vjust = 0.5, hjust = 1.1, size = 2, angle = 90)  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"), legend.position = "none", axis.text.y = element_text(angle = 0)) + scale_color_manual(values = c("darkorange3", "chartreuse4", "darkslategrey", "darkred"))  + coord_cartesian(ylim = c(-1, 1.5)) 
dev.off()

###################
# Figure 5
# Calculate Bray-Curtis Similarity for Trout Bog vs Mary Lake hypolimnia, 2007

# Select data
bog1 <- "TBH"
bog2 <- "MAH"
table <- clade_table07
query <- bog_subset(bog1, table)
database <- bog_subset(bog2, table)
# Remove January samples
query <- query[, c(1:32, 35:80)]

# Calculate Bray-Curtis Similarity between every dimictic and meromictic sample in the given subset, then average by dimictic sample
mean.bc <- c()
for (i in 1:dim(query)[2]) {
  bc.dis <- c()
  for (j in 1:dim(database)[2]) {
    test <- cbind(query[, i], database[, j])
    # Subtract from 1 to convert dissimilarity to similarity
    bc.dis[j] <- 1 - vegdist(t(test), method = "bray")
  }
  mean.bc[i] <- mean(bc.dis) 
}

# Make dataframe for plotting
dates <- extract_date(colnames(query))
plot.data <- data.frame(dates, mean.bc)
colnames(plot.data) <- c("Dates", "BrayCurtis")

TBHmat <- make_temp_matrix("TBH.....07", metadata)
TBHmat <- melt(TBHmat)
colnames(TBHmat) <- c("Depth", "Date", "Temperature")
TBHmat$Date <- as.Date(TBHmat$Date, format = "%Y-%m-%d")
TBHmat$Depth <- -TBHmat$Depth / 18.04 + 0.388

# Need to close polygons - add 0 or max values at top and bottom of graph
TBHmat <- TBHmat[which(is.na(TBHmat$Temperature) == F), ]
# For each date, add a hidden value outside of the plotting range
add <- TBHmat[which(TBHmat$Depth == 0.388), ]
add$Depth <- rep(-1, length(add$Depth))
add$Temperature <- rep(4, length(add$Temperature))
add2 <- add
add2$Depth <- rep(3, length(add2$Depth))
add2$Temperature <- rep(28, length(add2$Temperature))
TBHmat2 <- rbind(TBHmat, add, add2)


pdf(file = paste(path2repo, "TBH_v_MAH_bray_curtis.pdf", sep = ""), width = 3.3125, height = 2.5)
ggplot() + stat_contour(data = TBHmat2, aes(y = Depth, x = Date, z = Temperature, fill = ..level..), geom = "polygon") + scale_fill_gradientn(colours = c("dodgerblue", "cyan", "green", "yellow", "red"), "Temp", limits = c(4, 28)) + geom_line(data = plot.data, aes(x = Dates, y = BrayCurtis), size = 1.5) + labs(y = "Sorenson Similarity Index", x = NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "dodgerblue3"), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title.x = element_text(size = 15, vjust = 0.2), axis.title.y = element_text(size = 12, vjust = 1.6), axis.text.y = element_text(colour = "black", size = 10)) + coord_cartesian(xlim = extract_date(c("TBH07Jun07", "TBH11Nov07")), ylim = c(0.07, 0.4))
dev.off()

# Select data
clade_table08 <- year_subset("08", clade_table)
bog1 <- "NSH"
bog2 <- "MAH"
table <- clade_table08
query <- bog_subset(bog1, table)
database <- bog_subset(bog2, table)

# Calculate Bray-Curtis Similarity between every dimictic and meromictic sample in the given subset, then average by dimictic sample
mean.bc <- c()
for(i in 1:dim(query)[2]){
  bc.dis <- c()
  for(j in 1:dim(database)[2]){
    test <- cbind(query[, i], database[, j])
    #Subtract from 1 to convert dissimilarity to similarity
    bc.dis[j] <- 1 - vegdist(t(test), method="bray")
  }
  mean.bc[i] <- mean(bc.dis) 
}

# Make dataframe for plotting
dates <- extract_date(colnames(query))
plot.data <- data.frame(dates, mean.bc)
colnames(plot.data) <- c("Dates", "BrayCurtis")

NSHmat <- make_temp_matrix("NSH.....08", metadata)
# Remove columns with NAs
NSHmat <- NSHmat[,c(1:14, 16:30, 32:54)]
NSHmat <- melt(NSHmat)
colnames(NSHmat) <- c("Depth", "Date", "Temperature")
NSHmat$Date <- as.Date(NSHmat$Date, format = "%Y-%m-%d")
NSHmat$Depth <- -NSHmat$Depth / 10.71 + 0.43
# Need to close polygons - add 0 or max values at top and bottom of graph
NSHmat <- NSHmat[which(is.na(NSHmat$Temperature) == F), ]
# For each date, add a hidden -1 value
add <- NSHmat[which(NSHmat$Depth == 0.43),]
add$Depth <- rep(-1, length(add$Depth))
add$Temperature <- rep(4, length(add$Temperature))
add2 <- add
add2$Depth <- rep(3, length(add2$Depth))
add2$Temperature <- rep(28, length(add2$Temperature))
NSHmat2 <- rbind(NSHmat, add, add2)
# Remove dates with missing depth measurements

pdf(file = paste(path2repo, "NSH_v_MAH_bray_curtis.pdf", sep = ""), width = 3.3125, height = 2.5)
ggplot() + stat_contour(data = NSHmat2, aes(y = Depth, x = Date, z = Temperature, fill = ..level..), geom = "polygon") + coord_cartesian(xlim = extract_date(c("NSH25May08", "NSH01Nov08")), ylim = c(0.05, 0.42)) + scale_fill_gradientn(colours = c("dodgerblue", "cyan", "green", "yellow", "red"), "Temp", limits = c(4, 28)) + geom_line(data = plot.data, aes(x = Dates, y = BrayCurtis), size = 1.5) + labs(y = "Sorenson Similarity Index", x = NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill="dodgerblue3"), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title.x = element_text(size = 15, vjust = 0.2), axis.title.y = element_text(size = 12, vjust = 1.6), axis.text.y = element_text(colour = "black", size = 10)) 
dev.off()

#################
# Figure 6A
# Indicator analysis of epilimnia vs hypolimnia habitat preference

# Make a table that has OTUs grouped at every taxonomic level possible. This is so that groups at different levels will compete in the analyis.
# For example, I want to know if clade acI-B is a better or worse indicator than its phylum, Actinobacteria, so I run the analysis on both levels at once.

# Rename the OTUs with their full taxonomic assignment
named_otu_table <- otu_table
fullnames <- c()
for(i in 1:dim(taxonomy)[1]){
  fullnames[i] <- paste(taxonomy[i,], collapse = ";")
}
fullnames <- make.unique(fullnames)
rownames(named_otu_table) <- fullnames

# Re-run clade and phylum tables so that names are no longer shortened.
clade_table <- combine_otus("Clade", otu_table, taxonomy)
phylum_table <- combine_otus("Phylum", otu_table, taxonomy)
lineage_table <- combine_otus("Lineage", otu_table, taxonomy)
order_table <- combine_otus("Order", otu_table, taxonomy)
class_table <- combine_otus("Class", otu_table, taxonomy)

# Combine the tables at all taxonomic levels into one giant table
full_table <- rbind(named_otu_table, clade_table, lineage_table, order_table, class_table, phylum_table)

# Remove unclassified groups - I'm not interested in these for this analysis
classified <- grep("unclassified", rownames(full_table))
classified1 <- grep("__$", rownames(full_table))
parsed_table <- full_table[-c(classified, classified1),]

# Identify mixing dates
metadata$Sample_Name <- as.character(metadata$Sample_Name)
metalakes <- substr(metadata$Sample_Name, start=1, stop=3)
metayears <- substr(metadata$Sample_Name, start=9, stop=10)
metaTBH <- metadata[which(metalakes == "TBH"), c(1,2,4)]
metaTBH <- dcast(metaTBH, Sample_Name~Depth, fun.aggregate=mean)
TBHmixes <- metaTBH$Sample_Name[which(metaTBH$"0.5" - metaTBH$"7" < 1)]

metaNSH <- metadata[which(metalakes == "NSH"), c(1,2,4)]
metaNSH <- dcast(metaNSH, Sample_Name~Depth, fun.aggregate=mean)
NSHmixes <- metaNSH$Sample_Name[which(metaNSH$"0.5" - metaNSH$"4" < 1)]

metaSSH <- metadata[which(metalakes == "SSH"), c(1,2,4)]
metaSSH <- dcast(metaSSH, Sample_Name~Depth, fun.aggregate=mean)
SSHmixes <- metaSSH$Sample_Name[which(metaSSH$"0.5" - metaSSH$"8" < 1)]

metaCBH <- metadata[which(metalakes == "CBH"), c(1,2,4)]
metaCBH <- dcast(metaCBH, Sample_Name~Depth, fun.aggregate=mean)
CBHmixes <- metaCBH$Sample_Name[which(metaCBH$"0.5" - metaCBH$"2" < 1)]

metaWSH <- metadata[which(metalakes == "WSH"), c(1,2,4)]
metaWSH <- dcast(metaWSH, Sample_Name~Depth, fun.aggregate=mean)
WSHmixes <- metaWSH$Sample_Name[which(metaWSH$"0.5" - metaWSH$"4" < 1)]

metaFBH <- metadata[which(metalakes == "FBH"), c(1,2,4)]
metaFBH <- dcast(metaFBH, Sample_Name~Depth, fun.aggregate=mean)
FBHmixes <- metaFBH$Sample_Name[which(metaFBH$"0.5" - metaFBH$"1.5" < 1)]

mixes <- c(CBHmixes, FBHmixes, WSHmixes, TBHmixes, NSHmixes, SSHmixes)

# Use the full dataset for epi vs hypo  (minus polymictic lakes and mixing dates)
# Note: not using "NS." because there are some "NSU" (layer unknown) samples from that lake
input_table <- bog_subset("NSE|NSH|SS.|TB.|HK.|MA.", parsed_table)
mixing_dates <- match(mixes, substr(colnames(input_table), start=1, stop=10))
remove <- mixing_dates[which(is.na(mixing_dates) == F)]
input_table <- input_table[,-remove]

# Keep only groups with abundances in the top quantile (75th or higher)
threshold <- quantile(rowSums(input_table))[4]
input_table <- input_table[which(rowSums(input_table) >= threshold),]

# Format table for input into indicspecies analysis
input_table <- t(input_table)
input_table <- as.data.frame(input_table)
# Group by layer identifier
sampleids <- rownames(input_table)
layer <- substr(sampleids, start = 3, stop = 3)
layerid <- c("E", "H")

epi_v_hypo <- c()
for(i in 1:length(layerid)){
  epi_v_hypo[which(layer == layerid[i])] <- i
}

# Run multipatt(), specifing the the index of choice is group-normalized correlation
clade_by_layer <- multipatt(x = input_table, cluster = epi_v_hypo, func = "r.g", control = how(nperm = 9999))

# Sort results
results <- clade_by_layer$sign
epi <- results[which(results$index == 1), ]
epi <- epi[order(epi$stat, decreasing=T), ]

# Manually pick the top 10. Because groups from the same phylogeny are competing, choose the best indicator (by correlation coefficient) for an evolutionary branch.
# For example, if phylum Actinobacteria is a better indicator than acI-B, report only Actinobacteria
# But if acI_B is the better indicator, report both.
epi.indicators <- epi[c(1, 4, 6, 13, 14, 18, 21, 25, 26, 30), ]

# Add a column of abundance as % community for these indicators
epi_table <- bog_subset("..E", t(input_table))
hits <- match(rownames(epi.indicators), rownames(epi_table))
epi.indicators$abundance <- rowSums(epi_table[hits,]) / sum(rowSums(epi_table)) * 100

# Format for plotting
rownames(epi.indicators) <- substr(rownames(epi.indicators), start = 13, stop = 100)
rownames(epi.indicators) <- gsub("p__|c__|o__|\\[|\\]", "", rownames(epi.indicators))
epi.indicators$groups <- factor(rownames(epi.indicators), levels=rownames(epi.indicators))

pdf(file = paste(path2repo, "Epi_indicators.pdf", sep = ""), width = 3.3125 * 2, height = 4.625 / 2)
ggplot(data = epi.indicators, aes(x = groups, y = abundance, fill = stat)) + geom_bar(stat = "identity", colour = "black") + coord_flip() + labs(x = NULL, y = "% of Community") + theme(axis.text.x = element_text(angle = 0, size = 6, colour = "black"), axis.title = element_text(size = 10, vjust = -0.5), axis.text.y = element_text(colour="black", size = 10), legend.text = element_text(size = 6), axis.ticks=element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_gradient(low = "cyan", high = "darkblue")
dev.off()

# Repeat with hypolimnion results
hypo <- results[which(results$index == 2),]

hypo <- hypo[order(hypo$stat, decreasing=T),]
hypo.indicators <- hypo[c(1, 2, 7, 8, 9, 10, 14, 16, 18, 22), ]

hypo_table <- bog_subset("..H", t(input_table))
hits <- match(rownames(hypo.indicators), rownames(hypo_table))
hypo.indicators$abundance <- rowSums(hypo_table[hits,])/sum(rowSums(hypo_table)) * 100
rownames(hypo.indicators) <- substr(rownames(hypo.indicators), start=13, stop = 150)
rownames(hypo.indicators) <- gsub("p__|c__|o__|g__|f__|;g__;s__1||\\[|\\]", "", rownames(hypo.indicators))
hypo.indicators$groups <- factor(rownames(hypo.indicators), levels=rownames(hypo.indicators))

pdf(file = paste(path2repo, "Hypo_indicators.pdf", sep = ""), width = 3.3125 * 2 + 2, height = 4.625 / 2)
ggplot(data = hypo.indicators, aes(x = groups, y = abundance, fill = stat)) + geom_bar(stat = "identity", colour = "black") + coord_flip() + labs(x = NULL, y = "% of Community") + theme(axis.text.x = element_text(angle = 0, size = 6, colour = "black"), axis.title = element_text(size = 10, vjust = -0.5), axis.text.y = element_text(colour = "black", size = 10), legend.text = element_text(size = 6), axis.ticks = element_line(colour = "black")) + scale_y_continuous(expand = c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_gradient(low = "cyan", high = "darkblue")
dev.off()
