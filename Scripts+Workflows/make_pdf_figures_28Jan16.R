#Script to generate figures for NTL-MO manuscript
#Plots are saved as pdfs. For many figures, further processing is performed in Adobe Illustrator
###########

#Load packages
library(OTUtable)
library(vegan)
library(ggplot2)
library(reshape2)
library(grid)
library(indicspecies)

#Load data from OTUtable
data(otu_table)
data(taxonomy)
data(metadata)

#Generate tables at the clade and phylum level
tribe_table <- combine_otus("Tribe", otu_table, taxonomy)
clade_table <- combine_otus("Clade", otu_table, taxonomy)
lineage_table <- combine_otus("Lineage", otu_table, taxonomy)
order_table <- combine_otus("Order", otu_table, taxonomy)
class_table <- combine_otus("Class", otu_table, taxonomy)
phylum_table <- combine_otus("Phylum", otu_table, taxonomy)

#Reduce names to last known taxonmic information
clade_table <- reduce_names(clade_table)
phylum_table <- reduce_names(phylum_table)

clade_table07 <- year_subset("07", clade_table)
#Set the multiplot function (from user on Stack Overflow). Used to output multiple plots in one pdf
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
#Figure 1
#List lake categories
lakes <- c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA")

#Split epilimnion and hypolimnion into separate tables
epilimnia <- bog_subset("..E", otu_table)
hypolimnia <- bog_subset("..H", otu_table)

#Calculate observed richness
epi.chao1 <- apply(epilimnia, 2, obs_richness)
hypo.chao1 <- apply(hypolimnia, 2, obs_richness)

#Extract sampling location from sample names
epi.lakes <- substr(names(epi.chao1), start=1, stop=2)
hypo.lakes <- substr(names(hypo.chao1), start=1, stop=2)

#Make dataframe for plotting
epi.data <- data.frame(epi.lakes, epi.chao1)
hypo.data <- data.frame(hypo.lakes, hypo.chao1)
epi.data$epi.lakes <- ordered(epi.data$epi.lakes, levels = lakes)
hypo.data$hypo.lakes <- ordered(hypo.data$hypo.lakes, levels = lakes)

#1A
pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/epi_boxplot.pdf", width = 3.3125, height = 3)
ggplot(data=epi.data, aes(y=epi.chao1, x=epi.lakes, fill=epi.lakes)) + geom_boxplot() + labs(y="Observed Richness", x = NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks=element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust=0.5, vjust=0.1), axis.text.y = element_text(colour="black", size = 10), legend.position = "none") + scale_fill_brewer(palette="Set3")
dev.off()

#1B
pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/hypo_boxplot.pdf", width = 3.3125, height = 3)
ggplot(data=hypo.data, aes(y=hypo.chao1, x=hypo.lakes, fill=hypo.lakes)) + geom_boxplot() + labs(y="Observed Richness", x = NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks=element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title.y = element_text(size = 10, hjust=0.5, vjust=0.1), axis.text.y = element_text(colour="black", size = 10), legend.position="none") + scale_fill_brewer(palette="Set3")
dev.off()

#Check significance (not output as pdf, indicated as symbols in Illustrator)

epi.aov <- aov(epi.chao1 ~ epi.lakes, epi.data)
hypo.aov <- aov(hypo.chao1 ~ hypo.lakes, hypo.data)
TukeyHSD(epi.aov, "epi.lakes")
TukeyHSD(hypo.aov, "hypo.lakes")
###################
#Figure 2
#Trout Bog, 2007

#Identify mixing dates (less than 1 degree of temperature difference between 0.5 meters and maximum sampling depth)
metalakes <- substr(metadata$Sample_Name, start=1, stop=3)
metayears <- substr(metadata$Sample_Name, start=9, stop=10)
metaTBH <- metadata[which(metalakes == "TBH" & metayears == "07"), c(1,2,4)]
metaTBH <- dcast(metaTBH, Sample_Name~Depth, fun.aggregate=mean)
TBHmixes <- extract_date(metaTBH$Sample_Name[which(metaTBH$"0.5" - metaTBH$"7" < 1)])

#Make dataset of Trout Bog hypolimion samples from 2007
hypo <- bog_subset(paste("TBH", sep = ""), otu_table)
hypo <- year_subset("07", hypo)
#Calculate observed richness
hypo.rich <- apply(hypo, 2, obs_richness)
#Extract sampling date from sample names
hypo.date <- extract_date(colnames(hypo))
#Remove January samples - large gap distracts in plot, and winter samples are not considered in this study
hypo.rich <- hypo.rich[c(1:32, 35:80)]
hypo.date <- hypo.date[c(1:32, 35:80)]

#Make dataframe for plotting
TB_richness <- data.frame(hypo.date, hypo.rich)
colnames(TB_richness) <- c("date", "richness")

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/richness_over_time1.pdf", width = 3.3125, height = 2.3125)
ggplot() + geom_line(data=TB_richness, aes(x=date, y=richness), size=1) + labs(title = "Trout Bog", x = NULL, y = "Observed Richness") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks=element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title = element_text(size = 10, vjust=2), axis.text.y = element_text(colour="black", size=10), plot.title = element_text(size=12, vjust = 2), legend.position="none") + geom_point(data=TB_richness[match(TBHmixes, TB_richness$date),], aes(x=date, y=richness), size=2, colour = "red")
dev.off()


#Repeat with North Sparkling, 2008
metaNSH <- metadata[which(metalakes == "NSH" & metayears == "08"), c(1,2,4)]
metaNSH <- dcast(metaNSH, Sample_Name~Depth, fun.aggregate=mean)
NSHmixes <- extract_date(metaNSH$Sample_Name[which(metaNSH$"0.5" - metaNSH$"4" < 1)])

hypo <- bog_subset(paste("NSH", sep = ""), otu_table)
hypo <- year_subset("08", hypo)
hypo.rich <- apply(hypo, 2, obs_richness)
hypo.date <- extract_date(colnames(hypo))

NS_richness <- data.frame(hypo.date, hypo.rich)
colnames(NS_richness) <- c("date", "richness")

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/richness_over_time2.pdf", width = 3.3125, height = 2.3125)
ggplot() + geom_line(data=NS_richness, aes(x=date, y=richness), size=1.2) + labs(title = "North Sparkling Bog", x = NULL, y = "Observed Richness") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks=element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title = element_text(size = 10, vjust=2), axis.text.y = element_text(colour="black", size=10), plot.title = element_text(size=12, vjust = 2), legend.position = "none") + geom_point(data=NS_richness[match(NSHmixes, NS_richness$date),], aes(x=date, y=richness), size=2, colour = "red")
dev.off()

###################
#Figure 3A
#Calculate Bray-Curtis Similarity for Trout Bog vs Mary Lake hypolimnia, 2007

#Select data
bog1 <- "TBH"
bog2 <- "MAH"
table <- clade_table07
query <- bog_subset(bog1, table)
database <- bog_subset(bog2, table)
#Remove January samples
query <- query[,c(1:32, 35:80)]

#Calculate Bray-Curtis Similarity between every dimictic and meromictic sample in the given subset, then average by dimictic sample
mean.bc <- c()
for(i in 1:dim(query)[2]){
  bc.dis <- c()
  for(j in 1:dim(database)[2]){
    test <- cbind(query[,i], database[,j])
    #Subtract from 1 to convert dissimilarity to similarity
    bc.dis[j] <- 1 - vegdist(t(test), method="bray")
  }
  mean.bc[i] <- mean(bc.dis) 
}

#Make dataframe for plotting
dates <- extract_date(colnames(query))
plot.data <- data.frame(dates, mean.bc)
colnames(plot.data) <- c("Dates", "BrayCurtis")

TBHmat <- make_temp_matrix("TBH.....07", metadata)
TBHmat <- melt(TBHmat)
colnames(TBHmat) <- c("Depth", "Date", "Temperature")
TBHmat$Date <- as.Date(TBHmat$Date, format = "%Y-%m-%d")
TBHmat$Depth <- -TBHmat$Depth/18.04 + 0.388

#Need to close polygons - add 0 or max values at top and bottom of graph
TBHmat <- TBHmat[which(is.na(TBHmat$Temperature) == F),]
#For each date, add a hidden -1 value
add <- TBHmat[which(TBHmat$Depth == 0.388),]
add$Depth <- rep(-1, length(add$Depth))
add$Temperature <- rep(4, length(add$Temperature))
add2 <- add
add2$Depth <- rep(3, length(add2$Depth))
add2$Temperature <- rep(28, length(add2$Temperature))
TBHmat2 <- rbind(TBHmat, add, add2)


pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/TBH_v_MAH_bray_curtis.pdf", width = 3.3125, height = 2.5)
ggplot() + stat_contour(data=TBHmat2, aes(y=Depth, x=Date, z=Temperature, fill=..level..), geom="polygon") + scale_fill_gradientn(colours = c("dodgerblue", "cyan", "green", "yellow", "red"), "Temp", limits=c(4, 28)) + geom_line(data=plot.data, aes(x=Dates, y=BrayCurtis), size=1.5) + labs(y = "Bray-Curtis Similarity", x=NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill="dodgerblue3"), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title.x = element_text(size = 15, vjust=0.2), axis.title.y=element_text(size=15, vjust=1.6), axis.text.y = element_text(colour="black", size=10)) + coord_cartesian(xlim = extract_date(c("TBH07Jun07", "TBH11Nov07")), ylim=c(0.07, 0.372))
dev.off()

#3B
#Select data
clade_table08 <- year_subset("08", clade_table)
bog1 <- "NSH"
bog2 <- "MAH"
table <- clade_table08
query <- bog_subset(bog1, table)
database <- bog_subset(bog2, table)

#Calculate Bray-Curtis Similarity between every dimictic and meromictic sample in the given subset, then average by dimictic sample
mean.bc <- c()
for(i in 1:dim(query)[2]){
  bc.dis <- c()
  for(j in 1:dim(database)[2]){
    test <- cbind(query[,i], database[,j])
    #Subtract from 1 to convert dissimilarity to similarity
    bc.dis[j] <- 1 - vegdist(t(test), method="bray")
  }
  mean.bc[i] <- mean(bc.dis) 
}

#Make dataframe for plotting
dates <- extract_date(colnames(query))
plot.data <- data.frame(dates, mean.bc)
colnames(plot.data) <- c("Dates", "BrayCurtis")

NSHmat <- make_temp_matrix("NSH.....08", metadata)
#Remove columns with NAs
NSHmat <- NSHmat[,c(1:14, 16:30, 32:54)]
NSHmat <- melt(NSHmat)
colnames(NSHmat) <- c("Depth", "Date", "Temperature")
NSHmat$Date <- as.Date(NSHmat$Date, format = "%Y-%m-%d")
NSHmat$Depth <- -NSHmat$Depth/10.71 + 0.43
#Need to close polygons - add 0 or max values at top and bottom of graph
NSHmat <- NSHmat[which(is.na(NSHmat$Temperature) == F),]
#For each date, add a hidden -1 value
add <- NSHmat[which(NSHmat$Depth == 0.43),]
add$Depth <- rep(-1, length(add$Depth))
add$Temperature <- rep(4, length(add$Temperature))
add2 <- add
add2$Depth <- rep(3, length(add2$Depth))
add2$Temperature <- rep(28, length(add2$Temperature))
NSHmat2 <- rbind(NSHmat, add, add2)
#Remove dates with missing depth measurements

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/NSH_v_MAH_bray_curtis.pdf", width = 3.3125, height = 2.5)
ggplot() + stat_contour(data=NSHmat2, aes(y=Depth, x=Date, z=Temperature, fill=..level..), geom="polygon") + coord_cartesian(xlim = extract_date(c("NSH25May08", "NSH01Nov08")), ylim=c(0.05, 0.42)) + scale_fill_gradientn(colours = c("dodgerblue", "cyan", "green", "yellow", "red"), "Temp", limits=c(4, 28)) + geom_line(data=plot.data, aes(x=Dates, y=BrayCurtis), size=1.5) + labs(y = "Bray-Curtis Similarity", x=NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill="dodgerblue3"), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title.x = element_text(size = 15, vjust=0.2), axis.title.y=element_text(size=15, vjust=1.6), axis.text.y = element_text(colour="black", size=10)) 
dev.off()
#############
#Figure 4 - Network analysis
all.network <- read.table(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Network_analysis/allsamples_network_28Jan16.txt", header=T)
TBH <- bog_subset("TBH", otu_table)
NSH <- bog_subset("NSH", otu_table)
MAH <- bog_subset("MAH", otu_table)

TBH.edges <- table(c(as.character(all.network$index1), as.character(all.network$index2)))
TBH.nodes <- match(names(TBH.edges), rownames(TBH))
TBH.conn <- TBH[TBH.nodes,]
TBH.metric <- colSums(sweep(TBH.conn, 1, TBH.edges, "*"))
TBH.dates <- extract_date(colnames(TBH))

NSH.edges <- table(c(as.character(all.network$index1), as.character(all.network$index2)))
NSH.nodes <- match(names(NSH.edges), rownames(NSH))
NSH.conn <- NSH[NSH.nodes,]
NSH.metric <- colSums(sweep(NSH.conn, 1, NSH.edges, "*"))
NSH.dates <- extract_date(colnames(NSH))

MAH.edges <- table(c(as.character(all.network$index1), as.character(all.network$index2)))
MAH.nodes <- match(names(MAH.edges), rownames(MAH))
MAH.conn <- MAH[MAH.nodes,]
MAH.metric <- colSums(sweep(MAH.conn, 1, MAH.edges, "*"))
MAH.dates <- extract_date(colnames(MAH))

all.metric <- c(TBH.metric, NSH.metric, MAH.metric)
all.dates <- c(TBH.dates, NSH.dates, MAH.dates)
lakekey <- c(rep("TBH", length(TBH.metric)), rep("NSH", length(NSH.metric)), rep("MAH", length(MAH.metric)))
plot.conn <- data.frame(lakekey, all.dates, all.metric)
colnames(plot.conn) <- c("Lake", "Date", "Connectivity")
pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/connectivity2007.pdf", width = 3.3125*2, height = 2.3)
ggplot(data=plot.conn, aes(x=Date, y=Connectivity, colour=Lake)) + geom_line(size=1) + scale_y_log10() + coord_cartesian(xlim=extract_date(c("TBH01Jun07", "TBH16Nov07"))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title.x = element_text(size = 12, vjust=0.3), axis.title.y=element_text(size=12, vjust=1.3), axis.text.y = element_text(colour="black", size=10))
dev.off()

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/connectivity2008.pdf", width = 3.3125*2, height = 2.3)
ggplot(data=plot.conn, aes(x=Date, y=Connectivity, colour=Lake)) + geom_line(size=1) + scale_y_log10() + coord_cartesian(xlim=extract_date(c("TBH15May08", "TBH18Nov08"))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title.x = element_text(size = 12, vjust=0.3), axis.title.y=element_text(size=12, vjust=1.3), axis.text.y = element_text(colour="black", size=10))
dev.off()


#Use GLM

summary(glm(Connectivity ~ Lake + Date, family = "gaussian", data = plot.conn))

#Try looking at date in each lake individually
summary(glm(Connectivity ~ Date, family = "gaussian", data = plot.conn[which(plot.conn$Lake == "TBH"),]))
summary(glm(Connectivity ~ Date, family = "gaussian", data = plot.conn[which(plot.conn$Lake == "NSH"),]))
summary(glm(Connectivity ~ Date, family = "gaussian", data = plot.conn[which(plot.conn$Lake == "MAH"),]))
#Date is now significant in all three, although far less significant in Mary
#################
#Figure 5A
#Indicator analysis of epilimnia vs hypolimnia habitat preference

#Make a table that has OTUs grouped at every taxonomic level possible. This is so that groups at different levels will compete in the analyis.
#For example, I want to know if clade acI-B is a better or worse indicator than its phylum, Actinobacteria, so I run the analysis on both levels at once.

#Rename the OTUs with their full taxonomic assignment
named_otu_table <- otu_table
fullnames <- c()
for(i in 1:dim(taxonomy)[1]){
  fullnames[i] <- paste(taxonomy[i,], collapse = ";")
}
fullnames <- make.unique(fullnames)
rownames(named_otu_table) <- fullnames

#Re-run clade and phylum tables so that names are no longer shortened.
clade_table <- combine_otus("Clade", otu_table, taxonomy)
phylum_table <- combine_otus("Phylum", otu_table, taxonomy)

#Combine the tables at all taxonomic levels into one giant table
full_table <- rbind(named_otu_table, tribe_table, clade_table, lineage_table, order_table, class_table, phylum_table)

#Remove unclassified groups - I'm not interested in these for this analysis
classified <- grep("unclassified", rownames(full_table))
classified1 <- grep("__$", rownames(full_table))
parsed_table <- full_table[-c(classified, classified1),]

#Identify mixing dates
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

#
#Use the full dataset for epi vs hypo  (minus polymictic lakes and mixing dates)
#Note: not using "NS." because there are some "NSU" (layer unknown) samples from that lake
input_table <- bog_subset("NSE|NSH|SS.|TB.|HK.|MA.", parsed_table)
mixing_dates <- match(mixes, substr(colnames(input_table), start=1, stop=10))
remove <- mixing_dates[which(is.na(mixing_dates) == F)]
input_table <- input_table[,-remove]

#Keep only groups with abundances in the top quantile (75th or higher)
threshold <- quantile(rowSums(input_table))[4]
input_table <- input_table[which(rowSums(input_table) >= threshold),]

#Format table for input into indicspecies analysis
input_table <- t(input_table)
input_table <- as.data.frame(input_table)
#Group by layer identifier
sampleids <- rownames(input_table)
layer <- substr(sampleids, start=3, stop=3)
layerid <- c("E", "H")

epi_v_hypo <- c()
for(i in 1:length(layerid)){
  epi_v_hypo[which(layer == layerid[i])] <- i
}

#Run multipatt(), specifing the the index of choice is group-normalized correlation
clade_by_layer <- multipatt(x = input_table, cluster = epi_v_hypo, func = "r.g", control = how(nperm = 9999))

#Sort results
results <- clade_by_layer$sign
epi <- results[which(results$index == 1),]
epi <- epi[order(epi$stat, decreasing=T),]

#Manually pick the top 10. Because groups from the same phylogeny are competing, choose the best indicator (by correlation coefficient) for an evolutionary branch.
#For example, if phylum Actinobacteria is a better indicator than acI-B, report only Actinobacteria
#But if acI_B is the better indicator, report both.
epi.indicators <- epi[c(1, 2, 4, 6, 11, 16, 17, 20, 26, 39), ]

#Add a column of abundance as % community for these indicators
epi_table <- bog_subset("..E", t(input_table))
hits <- match(rownames(epi.indicators), rownames(epi_table))
epi.indicators$abundance <- rowSums(epi_table[hits,])/sum(rowSums(epi_table)) * 100

#Format for plotting
rownames(epi.indicators) <- substr(rownames(epi.indicators), start=13, stop = 100)
rownames(epi.indicators) <- gsub("p__|c__|o__|\\[|\\]", "", rownames(epi.indicators))
epi.indicators$groups <- factor(rownames(epi.indicators), levels=rownames(epi.indicators))

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/Epi_indicators.pdf", width = 3.3125*2, height = 4.625/2)
ggplot(data=epi.indicators, aes(x=groups, y=abundance, fill=stat)) + geom_bar(stat="identity", colour="black") + coord_flip() + labs(x = NULL, y = "% of Community") + theme(axis.text.x = element_text(angle = 0, size = 6, colour = "black"), axis.title = element_text(size = 10, vjust=-0.5), axis.text.y = element_text(colour="black", size = 6), legend.text = element_text(size=6), axis.ticks=element_line(colour="black")) + scale_y_continuous(expand = c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_gradient(low="cyan", high="darkblue")
dev.off()

#Repeat with hypolimnion results
hypo <- results[which(results$index == 2),]

hypo <- hypo[order(hypo$stat, decreasing=T),]
hypo.indicators <- hypo[c(1, 3, 4, 5, 6, 7, 10, 17, 21, 22), ]

hypo_table <- bog_subset("..H", t(input_table))
hits <- match(rownames(hypo.indicators), rownames(hypo_table))
hypo.indicators$abundance <- rowSums(hypo_table[hits,])/sum(rowSums(hypo_table)) * 100
rownames(hypo.indicators) <- substr(rownames(hypo.indicators), start=13, stop = 150)
rownames(hypo.indicators) <- gsub("p__|c__|o__|g__|f__|;g__;s__1||\\[|\\]", "", rownames(hypo.indicators))
hypo.indicators$groups <- factor(rownames(hypo.indicators), levels=rownames(hypo.indicators))

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/Hypo_indicators.pdf", width = 3.3125*2, height = 4.625/2)
ggplot(data=hypo.indicators, aes(x=groups, y=abundance, fill=stat)) + geom_bar(stat="identity", colour="black") + coord_flip() + labs(x = NULL, y = "% of Community") + theme(axis.text.x = element_text(angle = 0, size = 6, colour = "black"), axis.title = element_text(size = 10, vjust=-0.5), axis.text.y = element_text(colour="black", size = 6), legend.text = element_text(size=6), axis.ticks=element_line(colour="black")) + scale_y_continuous(expand = c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_fill_gradient(low="cyan", high="darkblue")
dev.off()
