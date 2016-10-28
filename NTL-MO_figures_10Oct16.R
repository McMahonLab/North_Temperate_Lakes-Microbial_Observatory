###Figures for NTL-MO paper, Oct 2016 edition###

#Set up the environment

library(OTUtable)
library(reshape2)
library(dplyr)
library(ggplot2)
library(indicspecies)
library(ape)
library(phyloseq)
library(vegan)
library(picante)
library(grid)
library(raster)
library(ggrepel)

data(metadata)

seq_table <- read.csv("C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_seqstable.csv", row.names = 1)
seq_taxonomy <- read.csv(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/seqs.98.cleantaxonomy.csv", row.names = 1, header = T, colClasses = c("character"))

seqs <- read.dna("C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/qc.bogs.clean.min25.fasta", format = "fasta")

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

#### Figure 1
# Read in fasta file of DNA sequences - needed for UniFrac distance calculation
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

# Set up phyloseq object to run UniFrac on
OTU <- otu_table(as.matrix(seq_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(seq_taxonomy))
sampledata <- sample_data(data.frame(Bog = substr(colnames(seq_table), start = 1, stop = 2), Layer = substr(colnames(seq_table), start = 3, stop = 3), Year = substr(colnames(seq_table), start = 9, stop = 10), row.names = colnames(seq_table), stringsAsfactors = F)) 
alldata <- phyloseq(OTU, TAX, sampledata, bogtree)

# Separate epilimnion and hypolimnion samples
epi <- prune_samples(sampledata$Layer == "E", alldata)
hypo <- prune_samples(sampledata$Layer == "H", alldata)

# Analyze and plot epilimnion points
x <- UniFrac(epi, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 1, stop = 3))
scores <- scores(pcoa)
years <- factor(substr(rownames(scores$sites), start = 9, stop = 10), levels = c("05", "07", "08", "09"))
lakes <- factor(substr(rownames(scores$sites), start = 1, stop = 2), levels = c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA"))
regime <- c()
regime[which(lakes == "CB" | lakes == "FB" | lakes == "WS")] <- "polymictic"
regime[which(lakes == "TB" | lakes == "NS" | lakes == "SS")] <- "dimictic"
regime[which(lakes == "HK" | lakes == "MA")] <- "meromictic"

plot.pcoa <- data.frame(scores$sites, years, lakes)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year", "Lake")

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

p1 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake, shape = Year)) + geom_point(size = 1, alpha = 1/2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), panel.border = element_rect(colour = "black", fill=NA, size=1)) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_shape_manual(values=c(21, 22, 23, 24))  + scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")

adonis(x ~ lakes, as(sample_data(epi), "data.frame"))
# r^2 = 0.34498, p = 0.001
# Calculate PERMADISP by mixing regime
adonis(x ~ regime, as(cbind(sample_data(epi), regime), "data.frame"))
# r^2 = 0.19574, p = 0.001

# Analyze and plot hypolimnion points
x <- UniFrac(hypo, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 1, stop = 3))
scores <- scores(pcoa)
years <- factor(substr(rownames(scores$sites), start = 9, stop = 10), levels = c("05", "07", "08", "09"))
lakes <- factor(substr(rownames(scores$sites), start = 1, stop = 2), levels = c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA"))
regime <- c()
regime[which(lakes == "CB" | lakes == "FB" | lakes == "WS")] <- "polymictic"
regime[which(lakes == "TB" | lakes == "NS" | lakes == "SS")] <- "polymictic"
regime[which(lakes == "HK" | lakes == "MA")] <- "meromictic"
plot.pcoa <- data.frame(scores$sites, years, lakes)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year", "Lake")
# Calculate percent variation explained of each axis
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

p2 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake, shape = Year)) + geom_point(size = 1, alpha = 1/2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), panel.border = element_rect(colour = "black", fill=NA, size=1)) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_shape_manual(values=c(21, 22, 23, 24)) + scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")

# Calculate significance level of clustering by group
adonis(x ~ lakes, as(sample_data(hypo), "data.frame"))
# r^2 = 0.48651, p = 0.001
# Calculate PERMADISP by mixing regime
adonis(x ~ regime, as(cbind(sample_data(hypo), regime), "data.frame"))
# r^2 = 0.21907, p = 0.001

#Reduce size of legends
pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figure1.pdf", width = 3.3125, height = 7)
multiplot(p1, p2, cols = 1)
dev.off()


###Figure 2
phylum_table <- combine_otus("Phylum", seq_table, seq_taxonomy)

abun <- phylum_table[which(rowSums(phylum_table) >= 10000),]
other <- phylum_table[which(rowSums(phylum_table) < 10000),]
new.phylum_table <- rbind(abun, colSums(other))

# Shorten up those names again and remove extraneous rownames
get.names <- strsplit(rownames(new.phylum_table), "p__")
phyla.names <- c()
for(i in 1:length(get.names)){
  phyla.names[i] <- get.names[[i]][2]
}
phyla.names[8] <- "unclassifed"
phyla.names[14] <- "other"
phyla.names <- factor(phyla.names, levels = rev(phyla.names))

# Convert to long format
new.phylum_table$seq_taxonomy <- phyla.names
new.phylum_table <- melt(new.phylum_table)

# Setup factors for plotting
new.phylum_table$Lake <- factor(substr(new.phylum_table$variable, start = 1, stop = 2), levels = c("CB", "FB", "NS", "WS", "TB", "SS", "HK", "MA"))
new.phylum_table$Layer <- substr(new.phylum_table$variable, start = 3, stop = 3)
new.phylum_table$LakeLayer <- factor(substr(new.phylum_table$variable, start = 1, stop = 3), levels = c("CBE", "FBE", "NSE", "WSE", "TBE", "SSE", "HKE", "MAE", "CBH", "FBH", "NSH", "WSH", "TBH", "SSH", "HKH", "MAH"))
new.phylum_table$Month <- factor(substr(new.phylum_table$variable, start = 6, stop = 8), levels = c("JUN", "JUL", "AUG", "SEP", "OCT","NOV", "JAN", "Feb", "MAR", "APR", "MAY"))
levels(new.phylum_table$Month)[levels(new.phylum_table$Month)=="Feb"] <- "FEB"

#levels(new.phylum_table$Month) <- c("JUN", "JUL", "AUG", "SEP", "OCT","NOV", "JAN", "FEB", "MAR", "APR", "MAY")
# Make color palette
pal2 = c("#005682", "#edfb48", "#a1a100", "#626262", "#008141", "#008282", "#00d5f2", "#f2a400", "#209401", "#929292", "#3885e7", "#ff8400", "#391826", "#f4bebe")

# Plot Trout Bog epi phyla
TBE_phyla <- filter(new.phylum_table, LakeLayer == "TBE")
TBE_phyla <- group_by(TBE_phyla, seq_taxonomy, Month)
TBE_phyla <- summarise(TBE_phyla, mean = mean(value))
p1 <- ggplot(data = TBE_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "Trout Bog Epilimnion") + theme(axis.text.x = element_text(size = 12, angle = 90, color="black"), axis.text.y = element_text(size=14, color="black"), axis.title = element_text(size = 15, vjust=2), legend.title = element_blank(), legend.text = element_text(size = 10)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

# Plot Trout Bog hypo phyla
TBH_phyla <- filter(new.phylum_table, LakeLayer == "TBH")
TBH_phyla <- group_by(TBH_phyla, seq_taxonomy, Month)
TBH_phyla <- summarise(TBH_phyla, mean = mean(value))
p2 <- ggplot(data = TBH_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "Trout Bog Hypolimnion") + theme(axis.text.x = element_text(size = 12, angle = 90, color="black"), axis.text.y = element_text(size=14, color="black"), axis.title = element_text(size = 15, vjust=2), legend.title = element_blank(), legend.text = element_text(size = 10)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figure2.pdf", width = 4.5, height = 7)
multiplot(p1, p2, cols=1)
dev.off()

### Figure 3
# Select lake and layer, then calculate UniFrac and plot vs time difference
TBH <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "H", alldata)
x <- UniFrac(TBH, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p1 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Trout Bog Hypolimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4))

TBE <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "E", alldata)
x <- UniFrac(TBE, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p2 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Trout Bog Epilimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4))

MAH <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer == "H", alldata)
x <- UniFrac(MAH, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p3 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Mary Lake Hypolimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4))

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figure3.pdf", width = 3.3125*2, height = 7)
multiplot(p1, p2, p3, cols = 1)
dev.off()

### Figure 4
# Pick lineages to investigate-in this case, PnecC and Rhodo in TBH
interest <- c("betIV$", "betI$")

# Make a plot of abundance of each OTU in the clades of interest over time in each year, then combine by each lineage at the end of the loop
for(i in 1:length(interest)){
  groups <- seq_table[grep(interest[i], seq_taxonomy$Lineage),]
  rownames(groups) <- make.unique(seq_taxonomy$Tribe[grep(interest[i], seq_taxonomy$Lineage)])
  TBH <- bog_subset("TBH", groups)
  TBH <- TBH[which(rowSums(TBH) > 1000), ]
  year1 <- year_subset("05", TBH)
  r <- cov2cor(cov(t(year1)))[2, 1]
  year1$seq <- rownames(year1)
  year1 <- melt(year1)
  year1$dates <- extract_date(year1$variable)
  p1 <- ggplot(data = year1, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + theme(legend.position = "none") + labs(title = "2005", x = "Date", y = "Abundance")  + xlim(extract_date(c("TBE15May05", "TBE20Nov05"))) + geom_label(data = NULL, x = as.numeric(extract_date(c("TBE01Oct05"))), y = max(year2$value) * 0.75, label = paste("r =", round(r, 2)), color = "black")

  
  
  year2 <- year_subset("07", TBH)
  r <- cov2cor(cov(t(year2)))[2, 1]
  year2$seq <- rownames(year2)
  year2 <- melt(year2)
  year2$dates <- extract_date(year2$variable)
  p2 <- ggplot(data = year2, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + theme(legend.position = "none") + labs(title = "2007", x = "Date", y = "Abundance") + xlim(extract_date(c("TBE15May07", "TBE20Nov07"))) + geom_label(data = NULL, x = as.numeric(extract_date(c("TBE01Oct07"))), y = max(year2$value) * 0.75, label = paste("r =", round(r, 2)), color = "black")
  
  year3 <- year_subset("08", TBH)
  r <- cov2cor(cov(t(year3)))[2, 1]
  year3$seq <- rownames(year3)
  year3 <- melt(year3)
  year3$dates <- extract_date(year3$variable)
  p3 <- ggplot(data = year3, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + theme(legend.position = "none") + labs(title = "2008", x = "Date", y = "Abundance")  + xlim(extract_date(c("TBE15May08", "TBE20Nov08"))) + geom_label(data = NULL, x = as.numeric(extract_date(c("TBE01Oct08"))), y = max(year3$value) * 0.75, label = paste("r =", round(r, 2)), color = "black")
  
  year4 <- year_subset("09", TBH)
  r <- cov2cor(cov(t(year4)))[2, 1]
  year4$seq <- rownames(year4)
  year4 <- melt(year4)
  year4$dates <- extract_date(year4$variable)
  p4 <- ggplot(data = year4, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + theme(legend.position = "none") +  labs(title = "2009", x = "Date", y = "Abundance")  + xlim(extract_date(c("TBE15May09", "TBE20Nov09"))) + geom_label(data = NULL, x = as.numeric(extract_date(c("TBE01Oct09"))), y = max(year4$value) * 0.75, label = paste("r =", round(r, 2)), color = "black")
  
  pdf(file = paste("C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figure4.", i, ".pdf", sep = ""), width = 5, height = 5)
  multiplot(p1, p2, p3, p4, cols = 2)
  dev.off()
}

### Figure 5
# Measure traits of lineages to show they are consistent
lineage_table <- combine_otus("Lineage", seq_table, seq_taxonomy)

interest <- c(";betII$", ";acI$", ";betI$", ";bacI$", ";gamI$", ";verI-A$", ";betIV$", ";acV$", ";gamIII$", ";alfI$", ";betIII$")

lineages <- lineage_table[grep(paste(interest, collapse = "|"), rownames(lineage_table)), ]
# I'm only showing a few lakes, but these can be adjusted easilty
lakes <- c("CBE", "TBE", "SSE", "MAE")

# For each lake, calculate my 3 metrics and plot. Save the plot to combine  into 1 pdf
for(i in 1:length(lakes)){
  table <- bog_subset(lakes[i], lineages)
  
  abundance <- c()
  persistence <- c()
  variance <- c()
  
  for(j in 1:dim(table)[1]){
    row <- table[j, ]
    abundance[j] <- sum(row[which(row > 0)])/length(which(row > 0))
    persistence[j] <- length(which(row > 0))/length(row)
    variance[j] <- cv(as.numeric(row))
  }
  
  to.plot <- data.frame(abundance, persistence, variance)
  to.plot$lineage <- rownames(reduce_names(table))
  
  z <- ggplot(data = to.plot, aes(x = variance, y = abundance/2500, fill = persistence)) + geom_point() + theme_bw() + labs(title = paste(lakes[i]), x = "Variability (CV)", y = "Mean Abundance when Present") + geom_label_repel(aes(label = lineage)) + scale_fill_gradient(low = "white", high = "lightgreen")

assign(paste("p", i, sep = ""), z)
}

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figure5.pdf", width = 3.3125*2.75, height = 7)
multiplot(p1, p2, p3, p4, cols = 2)
dev.off()

####Supplemental Figures####
### Figure S1 - richness boxplots
lakes <- c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA")

# Split epilimnion and hypolimnion into separate tables
epilimnia <- bog_subset("..E", seq_table)
hypolimnia <- bog_subset("..H", seq_table)

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
p1 <- ggplot(data = epi.data, aes(y = epi.obs, x = epi.lakes, fill = epi.lakes)) + geom_boxplot() + labs(y = "Observed Richness", x = NULL, title = "Epilimnion") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), legend.position = "none") + scale_fill_brewer(palette = "Set3") + theme_bw()

# plot 1B
p2 <- ggplot(data = hypo.data, aes(y = hypo.obs, x = hypo.lakes, fill = hypo.lakes)) + geom_boxplot() + labs(y="Observed Richness", x = NULL, title = "Hypolimnion") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title.y = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), legend.position="none") + scale_fill_brewer(palette = "Set3") + theme_bw()

# Check significance using the Wilcoxon Rank Test on medians (not output as pdf, indicated as symbols in Illustrator)
pairwise.wilcox.test(epi.data$epi.obs, epi.data$epi.lakes, p.adjust.method = "bonferroni")

pairwise.wilcox.test(hypo.data$hypo.obs, hypo.data$hypo.lakes, p.adjust.method = "bonferroni")

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS1.pdf", width = 3.3125*2, height = 7)
multiplot(p1, p2, cols = 1)
dev.off()

### Figure S2 - richness over time
# Trout Bog, 2007
# Identify mixing dates (less than 1 degree of temperature difference between 0.5 meters and maximum sampling depth)
metalakes <- substr(metadata$Sample_Name, start = 1, stop = 3)
metayears <- substr(metadata$Sample_Name, start = 9, stop = 10)
metaTBH <- metadata[which(metalakes == "TBH" & metayears == "07"), c(1,2,4)]
metaTBH <- dcast(metaTBH, Sample_Name ~ Depth, fun.aggregate = mean)
TBHmixes <- extract_date(metaTBH$Sample_Name[which(metaTBH$"0.5" - metaTBH$"7" < 1)])

# Make dataset of Trout Bog hypolimion samples from 2007
hypo <- bog_subset(paste("TBH", sep = ""), seq_table)
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

p1 <- ggplot() + geom_line(data = TB_richness, aes(x = date, y = richness), size = 1) + labs(title = "Trout Bog", x = NULL, y = "Observed Richness") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks=element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title = element_text(size = 10, vjust = 2), axis.text.y = element_text(colour = "black", size=10), plot.title = element_text(size = 12, vjust = 2), legend.position = "none") + geom_point(data = TB_richness[match(TBHmixes, TB_richness$date), ], aes(x = date, y = richness), size = 2, colour = "red") + theme_bw()
dev.off()


# Repeat with North Sparkling, 2008
metaNSH <- metadata[which(metalakes == "NSH" & metayears == "08"), c(1,2,4)]
metaNSH <- dcast(metaNSH, Sample_Name ~ Depth, fun.aggregate = mean)
NSHmixes <- extract_date(metaNSH$Sample_Name[which(metaNSH$"0.5" - metaNSH$"4" < 1)])

hypo <- bog_subset(paste("NSH", sep = ""), seq_table)
hypo <- year_subset("08", hypo)
hypo.rich <- apply(hypo, 2, obs_richness)
hypo.date <- extract_date(colnames(hypo))

NS_richness <- data.frame(hypo.date, hypo.rich)
colnames(NS_richness) <- c("date", "richness")

p2 <- ggplot() + geom_line(data = NS_richness, aes(x = date, y = richness), size = 1.2) + labs(title = "North Sparkling Bog", x = NULL, y = "Observed Richness") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title = element_text(size = 10, vjust=2), axis.text.y = element_text(colour = "black", size = 10), plot.title = element_text(size=12, vjust = 2), legend.position = "none") + geom_point(data = NS_richness[match(NSHmixes, NS_richness$date), ], aes(x = date, y = richness), size = 2, colour = "red") + theme_bw()

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS2.pdf", width = 3.3125*2, height = 5.5)
multiplot(p1, p2, cols = 1)
dev.off()

# Identify mixing dates (less than 1 degree of temperature difference between 0.5 meters and maximum sampling depth)
metalakes <- substr(metadata$Sample_Name, start = 1, stop = 3)
metayears <- substr(metadata$Sample_Name, start = 9, stop = 10)
metaTBH <- metadata[which(metalakes == "TBH" & metayears == "07"), c(1,2,4)]
metaTBH <- dcast(metaTBH, Sample_Name~Depth, fun.aggregate=mean)
TBHmixes <- extract_date(metaTBH$Sample_Name[which(metaTBH$"0.5" - metaTBH$"7" < 1)])

# Make dataset of Trout Bog hypolimion samples from 2007
hypo <- bog_subset(paste("TBH", sep = ""), seq_table)
hypo <- year_subset("07", hypo)
#Calculate observed richness
hypo.even <- apply(hypo, 2, pielou)
#Extract sampling date from sample names
hypo.date <- extract_date(colnames(hypo))
#Remove January samples - large gap distracts in plot, and winter samples are not considered in this study
hypo.even <- hypo.even[c(1:32, 35:80)]
hypo.date <- hypo.date[c(1:32, 35:80)]

# Make dataframe for plotting
TB_evenness <- data.frame(hypo.date, hypo.even)
colnames(TB_evenness) <- c("date", "evenness")

p1 <- ggplot() + geom_line(data = TB_evenness, aes(x = date, y = evenness), size = 1) + labs(title = "Trout 2007", x = NULL, y = "Pielou's Evenness") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks=element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title = element_text(size = 10, vjust = 2), axis.text.y = element_text(colour = "black", size=10), plot.title = element_text(size = 12, vjust = 2), legend.position = "none") + geom_point(data = TB_evenness[match(TBHmixes, TB_evenness$date), ], aes(x = date, y = evenness), size = 2, colour = "red") + theme_bw()

metalakes <- substr(metadata$Sample_Name, start = 1, stop = 3)
metayears <- substr(metadata$Sample_Name, start = 9, stop = 10)
metaNSH <- metadata[which(metalakes == "NSH" & metayears == "08"), c(1,2,4)]
metaNSH <- dcast(metaNSH, Sample_Name~Depth, fun.aggregate=mean)
NSHmixes <- extract_date(metaNSH$Sample_Name[which(metaNSH$"0.5" - metaNSH$"4" < 1)])

# Make dataset of Trout Bog hypolimion samples from 2007
hypo <- bog_subset(paste("NSH", sep = ""), seq_table)
hypo <- year_subset("08", hypo)
#Calculate observed richness
hypo.even <- apply(hypo, 2, pielou)
#Extract sampling date from sample names
hypo.date <- extract_date(colnames(hypo))
#Remove January samples - large gap distracts in plot, and winter samples are not considered in this study
hypo.even <- hypo.even[c(1:32, 35:80)]
hypo.date <- hypo.date[c(1:32, 35:80)]

# Make dataframe for plotting
NS_evenness <- data.frame(hypo.date, hypo.even)
colnames(NS_evenness) <- c("date", "evenness")

p2 <- ggplot() + geom_line(data = NS_evenness, aes(x = date, y = evenness), size = 1) + labs(title = "North Sparkling 2008", x = NULL, y = "Pielou's Evenness") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks=element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.title = element_text(size = 10, vjust = 2), axis.text.y = element_text(colour = "black", size=10), plot.title = element_text(size = 12, vjust = 2), legend.position = "none") + geom_point(data = NS_evenness[match(NSHmixes, NS_evenness$date), ], aes(x = date, y = evenness), size = 2, colour = "red") + theme_bw()


pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS2.1.pdf", width = 3.3125*2, height = 5.5)
multiplot(p1, p2, cols = 1)
dev.off()

### Figure S3 - Phylum comp by lake

# Specify the sampled layers
layer <- c("CBE", "FBE", "WSE", "NSE", "TBE", "SSE", "HKE", "MAE", "CBH", "FBH", "WSH", "NSH", "TBH",  "SSH", "HKH", "MAH")

# Set up a dataframe with the totals of each phylum for each sampling site
phylum_table <- combine_otus("Phylum", seq_table, seq_taxonomy)
layer.phyla <- rep(NA, dim(phylum_table)[1])

for(i in 1:length(layer)){
  dataset <- bog_subset(layer[i], phylum_table)
  layer.phyla <- cbind(layer.phyla, rowSums(dataset)) 
}

layer.phyla <- layer.phyla[,2:dim(layer.phyla)[2]]
colnames(layer.phyla) <- layer

# Combine low abundance phyla into a single category called "other"
abun <- layer.phyla[which(rowSums(layer.phyla) >= 10000),]
other <- layer.phyla[which(rowSums(layer.phyla) < 10000),]
new.layer <- rbind(abun, colSums(other))

# Shorten up those names again and remove extraneous rownames
get.names <- strsplit(rownames(new.layer), "p__")
phyla.names <- c()
for(i in 1:length(get.names)){
  phyla.names[i] <- get.names[[i]][2]
}
phyla.names[8] <- "unclassified"
phyla.names[14] <- "other"
phyla.names <- factor(phyla.names, levels = rev(phyla.names))
rownames(new.layer) <- NULL

# Convert data into a long format dataframe for use in ggplot
phyla_by_bog <- data.frame(phyla.names, new.layer)
phyla_by_bog2 <- melt(phyla_by_bog)

# Create color palette that can handle the large number of categories
pal2 = c("#005682", "#edfb48", "#a1a100", "#626262", "#008141", "#008282", "#00d5f2", "#f2a400", "#209401", "#929292", "#3885e7", "#ff8400", "#391826", "#f4bebe")

# Plot data as a stacked barplot
pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS3.pdf", width = 3.3125*2, height = 6)
ggplot(data=phyla_by_bog2, aes(x=variable, y=value, fill=phyla.names)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads") + theme(axis.text.x = element_text(size = 12, angle = 90, color="black"), axis.text.y = element_text(size=14, color="black"), axis.title = element_text(size = 15, vjust=2), legend.title = element_blank(), legend.text = element_text(size = 16)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 
dev.off()

### Figure S4 - Phylum comp by season by more lakes?

phylum_table <- combine_otus("Phylum", seq_table, seq_taxonomy)

abun <- phylum_table[which(rowSums(phylum_table) >= 10000),]
other <- phylum_table[which(rowSums(phylum_table) < 10000),]
new.phylum_table <- rbind(abun, colSums(other))

# Shorten up those names again and remove extraneous rownames
get.names <- strsplit(rownames(new.phylum_table), "p__")
phyla.names <- c()
for(i in 1:length(get.names)){
  phyla.names[i] <- get.names[[i]][2]
}
phyla.names[8] <- "unclassifed"
phyla.names[14] <- "other"
phyla.names <- factor(phyla.names, levels = rev(phyla.names))

# Convert to long format
new.phylum_table$seq_taxonomy <- phyla.names
new.phylum_table <- melt(new.phylum_table)

# Setup factors for plotting
new.phylum_table$Lake <- factor(substr(new.phylum_table$variable, start = 1, stop = 2), levels = c("CB", "FB", "NS", "WS", "TB", "SS", "HK", "MA"))
new.phylum_table$Layer <- substr(new.phylum_table$variable, start = 3, stop = 3)
new.phylum_table$LakeLayer <- factor(substr(new.phylum_table$variable, start = 1, stop = 3), levels = c("CBE", "FBE", "NSE", "WSE", "TBE", "SSE", "HKE", "MAE", "CBH", "FBH", "NSH", "WSH", "TBH", "SSH", "HKH", "MAH"))
new.phylum_table$Month <- factor(substr(new.phylum_table$variable, start = 6, stop = 8), levels = c("JAN", "Feb", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT","NOV"))
levels(new.phylum_table$Month) <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT","NOV")

# Make color palette
pal2 = c("#005682", "#edfb48", "#a1a100", "#626262", "#008141", "#008282", "#00d5f2", "#f2a400", "#209401", "#929292", "#3885e7", "#ff8400", "#391826", "#f4bebe")

CBE_phyla <- filter(new.phylum_table, LakeLayer == "CBE")
CBE_phyla <- group_by(CBE_phyla, seq_taxonomy, Month)
CBE_phyla <- summarise(CBE_phyla, mean = mean(value))
p1 <- ggplot(data = CBE_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "Crystal Bog Epilimnion") + theme(axis.text.x = element_text(size = 10, angle = 90, color="black"), axis.text.y = element_text(size=12, color="black"), axis.title = element_text(size = 12, vjust=2), legend.position = "none", plot.title = element_text(size = 12)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

CBH_phyla <- filter(new.phylum_table, LakeLayer == "CBH")
CBH_phyla <- group_by(CBH_phyla, seq_taxonomy, Month)
CBH_phyla <- summarise(CBH_phyla, mean = mean(value))
p2 <- ggplot(data = CBH_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "Crystal Bog Hypolimnion") + theme(axis.text.x = element_text(size = 10, angle = 90, color="black"), axis.text.y = element_text(size=12, color="black"), axis.title = element_text(size = 12, vjust=2) , legend.position = "none", plot.title = element_text(size = 12)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

NSE_phyla <- filter(new.phylum_table, LakeLayer == "NSE")
NSE_phyla <- group_by(NSE_phyla, seq_taxonomy, Month)
NSE_phyla <- summarise(NSE_phyla, mean = mean(value))
p3 <- ggplot(data = NSE_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "North Sparkling Epilimnion") + theme(axis.text.x = element_text(size = 10, angle = 90, color="black"), axis.text.y = element_text(size=12, color="black"), axis.title = element_text(size = 12, vjust=2), legend.position = "none", plot.title = element_text(size = 12)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

NSH_phyla <- filter(new.phylum_table, LakeLayer == "NSH")
NSH_phyla <- group_by(NSH_phyla, seq_taxonomy, Month)
NSH_phyla <- summarise(NSH_phyla, mean = mean(value))
p4 <- ggplot(data = NSH_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "North Sparkling Hypolimnion") + theme(axis.text.x = element_text(size = 10, angle = 90, color="black"), axis.text.y = element_text(size=10, color="black"), axis.title = element_text(size = 12, vjust=2), legend.position = "none", plot.title = element_text(size = 12)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS4.1.pdf", width = 5.5, height = 7)
multiplot(p1, p2, p3, p4, cols=2)
dev.off()

# Panel 2
SSE_phyla <- filter(new.phylum_table, LakeLayer == "SSE")
SSE_phyla <- group_by(SSE_phyla, seq_taxonomy, Month)
SSE_phyla <- summarise(SSE_phyla, mean = mean(value))
p1 <- ggplot(data = SSE_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "South Sparkling Epilimnion") + theme(axis.text.x = element_text(size = 10, angle = 90, color="black"), axis.text.y = element_text(size=12, color="black"), axis.title = element_text(size = 12, vjust=2), legend.position = "none", plot.title = element_text(size = 12)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

SSH_phyla <- filter(new.phylum_table, LakeLayer == "SSH")
SSH_phyla <- group_by(SSH_phyla, seq_taxonomy, Month)
SSH_phyla <- summarise(SSH_phyla, mean = mean(value))
p2 <- ggplot(data = SSH_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "South Sparkling Hypolimnion") + theme(axis.text.x = element_text(size = 10, angle = 90, color="black"), axis.text.y = element_text(size=12, color="black"), axis.title = element_text(size = 12, vjust=2) , legend.position = "none", plot.title = element_text(size = 12)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

MAE_phyla <- filter(new.phylum_table, LakeLayer == "MAE")
MAE_phyla <- group_by(MAE_phyla, seq_taxonomy, Month)
MAE_phyla <- summarise(MAE_phyla, mean = mean(value))
p3 <- ggplot(data = MAE_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "Mary Lake Epilimnion") + theme(axis.text.x = element_text(size = 10, angle = 90, color="black"), axis.text.y = element_text(size=12, color="black"), axis.title = element_text(size = 12, vjust=2), legend.position = "none", plot.title = element_text(size = 12)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

MAH_phyla <- filter(new.phylum_table, LakeLayer == "MAH")
MAH_phyla <- group_by(MAH_phyla, seq_taxonomy, Month)
MAH_phyla <- summarise(MAH_phyla, mean = mean(value))
p4 <- ggplot(data = MAH_phyla, aes(x = Month, y = mean, fill = seq_taxonomy)) + geom_bar(stat="identity", position = "fill") + labs(x = NULL, y = "Proportion of Observed Reads", title = "Mary Lake Hypolimnion") + theme(axis.text.x = element_text(size = 10, angle = 90, color="black"), axis.text.y = element_text(size=10, color="black"), axis.title = element_text(size = 12, vjust=2), legend.position = "none", plot.title = element_text(size = 12)) + scale_fill_manual(values=rev(pal2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"))  + scale_y_continuous(expand = c(0,0)) 

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS4.2.pdf", width = 5.5, height = 7)
multiplot(p1, p2, p3, p4, cols=2)
dev.off()


### Figure S5 - MNTD lengths

# I already measured MNTD for each sample using the following code and saved it to a file:
# mntd <- ses.mntd(t(seq_table), d, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
# write.csv(mntd, file = "mntd_taxon_relatedness_index.csv")

#Let's read it back in and start from there
mntd <- read.csv("C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/mntd_taxon_relatedness_index.csv", row.names = 1)
mntd$nti <- (mntd$mntd.rand.mean - mntd$mntd.obs)/mntd$mntd.rand.sd

mntd.results <- data.frame(mntd$mntd.obs.z)
mntd.results$site <- substr(rownames(mntd), start = 1, stop = 3)
mntd.results$site <- factor(mntd.results$site, levels = c("CBE", "FBE", "WSE", "NSE", "TBE", "SSE", "HKE", "MAE", "CBH", "FBH", "WSH", "NSH", "TBH", "SSH", "HKH", "MAH"))

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS5.pdf", width = 3.3125*2, height = 4)
ggplot(data = mntd.results[which(is.na(mntd.results$site) == F), ], aes(y = mntd.mntd.obs.z, x = site)) + geom_boxplot() + theme_bw() + labs(title = "Phylogenetic Clustering", x = NULL, y = "Mean Nearest Taxon Distance")
dev.off()

pairwise.wilcox.test(mntd.results$mntd.mntd.obs.z, mntd.results$site, p.adjust.method = "bonferroni")

### Figure S6 - More examples of OTU trends
interest <- c("^acI$")

for(i in 1:length(interest)){
  groups <- seq_table[grep(interest[i], seq_taxonomy$Lineage),]
  rownames(groups) <- make.unique(seq_taxonomy$Tribe[grep(interest[i], seq_taxonomy$Lineage)])
  TBE <- bog_subset("TBE", groups)
  TBE <- TBE[which(rowSums(TBE) > 1000), ]
  if(dim(TBE)[1] > 1){
    year1 <- year_subset("05", TBE)
    # print(cov2cor(cov(t(year1))))
    year1$seq <- rownames(year1)
    year1 <- melt(year1)
    year1$dates <- extract_date(year1$variable)
    p1 <- ggplot(data = year1, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2005 TBE", x = "Date", y = "Abundance")  + xlim(extract_date(c("TBE15May05", "TBE20Nov05"))) + theme(legend.position = "none")
    
    
    year2 <- year_subset("07", TBE)
    # print(cov2cor(cov(t(year2))))
    year2$seq <- rownames(year2)
    year2 <- melt(year2)
    year2$dates <- extract_date(year2$variable)
    p2 <- ggplot(data = year2, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2007", x = "Date", y = "Abundance") + xlim(extract_date(c("TBE15May07", "TBE20Nov07")))  + theme(legend.position = "none")
    
    year3 <- year_subset("08", TBE)
    # print(cov2cor(cov(t(year3))))
    year3$seq <- rownames(year3)
    year3 <- melt(year3)
    year3$dates <- extract_date(year3$variable)
    p3 <- ggplot(data = year3, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2008", x = "Date", y = "Abundance")  + xlim(extract_date(c("TBE15May08", "TBE20Nov08"))) + theme(legend.position = "none")  
    
    year4 <- year_subset("09", TBE)
    # print(cov2cor(cov(t(year4))))
    year4$seq <- rownames(year4)
    year4 <- melt(year4)
    year4$dates <- extract_date(year4$variable)
    p4 <- ggplot(data = year4, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2009", x = "Date", y = "Abundance")  + xlim(extract_date(c("TBE15May09", "TBE20Nov09"))) + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title = element_blank())
    
  
    pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS6.1.pdf", width = 3.3125*2, height = 6)
    multiplot(p1, p2, p3, p4, cols = 2)
    dev.off()
    
  }
  
}

#Try with MAE

interest <- c("^acI$", "betII$")

for(i in 1:length(interest)){
  groups <- seq_table[grep(interest[i], seq_taxonomy$Lineage),]
  rownames(groups) <- make.unique(seq_taxonomy$Tribe[grep(interest[i], seq_taxonomy$Lineage)])
  MAE <- bog_subset("MAE", groups)
  MAE <- MAE[which(rowSums(MAE) > 1000), ]
  if(dim(MAE)[1] > 1){
    year1 <- year_subset("05", MAE)
    # print(cov2cor(cov(t(year1))))
    year1$seq <- rownames(year1)
    year1 <- melt(year1)
    year1$dates <- extract_date(year1$variable)
    p1 <- ggplot(data = year1, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2005 MAE", x = "Date", y = "Abundance")  + xlim(extract_date(c("MAE15May05", "MAE20Nov05")))  + theme(legend.position = "none")  
    
    
    year2 <- year_subset("07", MAE)
    # print(cov2cor(cov(t(year2))))
    year2$seq <- rownames(year2)
    year2 <- melt(year2)
    year2$dates <- extract_date(year2$variable)
    p2 <- ggplot(data = year2, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2007", x = "Date", y = "Abundance") + xlim(extract_date(c("MAE15May07", "MAE20Nov07")))  + theme(legend.position = "none")  
    
    year3 <- year_subset("08", MAE)
    # print(cov2cor(cov(t(year3))))
    year3$seq <- rownames(year3)
    year3 <- melt(year3)
    year3$dates <- extract_date(year3$variable)
    p3 <- ggplot(data = year3, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2008", x = "Date", y = "Abundance")  + xlim(extract_date(c("MAE15May08", "MAE20Nov08")))  + theme(legend.position = "none")  
    
    year4 <- year_subset("09", MAE)
    # print(cov2cor(cov(t(year4))))
    year4$seq <- rownames(year4)
    year4 <- melt(year4)
    year4$dates <- extract_date(year4$variable)
    p4 <- ggplot(data = year4, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2009", x = "Date", y = "Abundance")  + xlim(extract_date(c("MAE15May09", "MAE20Nov09")))  + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title = element_blank())
    
    pdf(file = paste("C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS6.2.", i, ".pdf", sep = ""), width = 3.3125*2, height = 6)
    multiplot(p1, p2, p3, p4, cols = 2)
    dev.off()
    
  }
  
}

interest <- c("betIV$", "betI$")

for(i in 1:length(interest)){
  groups <- seq_table[grep(interest[i], seq_taxonomy$Lineage),]
  rownames(groups) <- make.unique(seq_taxonomy$Tribe[grep(interest[i], seq_taxonomy$Lineage)])
  TBH <- bog_subset("TBH", groups)
  TBH <- TBH[which(rowSums(TBH) > 1000), ]
  year1 <- year_subset("05", TBH)
  # print(cov2cor(cov(t(year1))))
  year1$seq <- rownames(year1)
  year1 <- melt(year1)
  year1$dates <- extract_date(year1$variable)
  p1 <- ggplot(data = year1, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2005 TBH", x = "Date", y = "Abundance")  + xlim(extract_date(c("TBE15May05", "TBE20Nov05"))) + theme(legend.position = "none")
  
  
  year2 <- year_subset("07", TBH)
  # print(cov2cor(cov(t(year2))))
  year2$seq <- rownames(year2)
  year2 <- melt(year2)
  year2$dates <- extract_date(year2$variable)
  p2 <- ggplot(data = year2, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2007", x = "Date", y = "Abundance") + xlim(extract_date(c("TBE15May07", "TBE20Nov07"))) + theme(legend.position = "none")
  
  year3 <- year_subset("08", TBH)
  # print(cov2cor(cov(t(year3))))
  year3$seq <- rownames(year3)
  year3 <- melt(year3)
  year3$dates <- extract_date(year3$variable)
  p3 <- ggplot(data = year3, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2008", x = "Date", y = "Abundance")  + xlim(extract_date(c("TBE15May08", "TBE20Nov08"))) + theme(legend.position = "none")
  
  year4 <- year_subset("09", TBH)
  # print(cov2cor(cov(t(year4))))
  year4$seq <- rownames(year4)
  year4 <- melt(year4)
  year4$dates <- extract_date(year4$variable)
  p4 <- ggplot(data = year4, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2009", x = "Date", y = "Abundance")  + xlim(extract_date(c("TBE15May09", "TBE20Nov09")))  + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title = element_blank())
  
  pdf(file = paste("C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS6.3.", i, ".pdf", sep = ""), width = 3.3125*2, height = 6)
  multiplot(p1, p2, p3, p4, cols = 2)
  dev.off()
  
}

#Try MAH

interest <- c("gamI$")

for(i in 1:length(interest)){
  groups <- seq_table[grep(interest[i], seq_taxonomy$Lineage),]
  rownames(groups) <- make.unique(seq_taxonomy$Tribe[grep(interest[i], seq_taxonomy$Lineage)])
  MAH <- bog_subset("MAH", groups)
  MAH <- MAH[which(rowSums(MAH) > 1000), ]
  if(dim(MAH)[1] > 1){
    year1 <- year_subset("05", MAH)
    # print(cov2cor(cov(t(year1))))
    year1$seq <- rownames(year1)
    year1 <- melt(year1)
    year1$dates <- extract_date(year1$variable)
    p1 <- ggplot(data = year1, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2005 MAH", x = "Date", y = "Abundance")  + xlim(extract_date(c("MAH15May05", "MAH20Nov05"))) + theme(legend.position = "none")
    
    
    year2 <- year_subset("07", MAH)
    # print(cov2cor(cov(t(year2))))
    year2$seq <- rownames(year2)
    year2 <- melt(year2)
    year2$dates <- extract_date(year2$variable)
    p2 <- ggplot(data = year2, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2007", x = "Date", y = "Abundance") + xlim(extract_date(c("MAH15May07", "MAH20Nov07"))) + theme(legend.position = "none")
    
    year3 <- year_subset("08", MAH)
    # print(cov2cor(cov(t(year3))))
    year3$seq <- rownames(year3)
    year3 <- melt(year3)
    year3$dates <- extract_date(year3$variable)
    p3 <- ggplot(data = year3, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2008", x = "Date", y = "Abundance")  + xlim(extract_date(c("MAH15May08", "MAH20Nov08"))) + theme(legend.position = "none")
    
    year4 <- year_subset("09", MAH)
    # print(cov2cor(cov(t(year4))))
    year4$seq <- rownames(year4)
    year4 <- melt(year4)
    year4$dates <- extract_date(year4$variable)
    p4 <- ggplot(data = year4, aes(x = dates, y = value, color = seq)) + geom_point() + geom_line() + theme_bw() + labs(title = "2009", x = "Date", y = "Abundance")  + xlim(extract_date(c("MAH15May09", "MAH20Nov09")))  + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.title = element_blank())
    
    pdf(file = paste("C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS6.4.", i, ".pdf", sep = ""), width = 3.3125*2, height = 6)
    multiplot(p1, p2, p3, p4, cols = 2)
    dev.off()
    
  }
  
}

### Figure S7 - Consistent lineage traits by year
lineage_table <- combine_otus("Lineage", seq_table, seq_taxonomy)

interest <- c(";betII$", ";acI$", ";betI$", ";bacI$", ";gamI$", ";verI-A$", ";betIV$", ";acV$", ";gamIII$", ";alfI$", ";betIII$")

lineages <- lineage_table[grep(paste(interest, collapse = "|"), rownames(lineage_table)), ]

year <- c("05", "07", "08", "09")
for(i in 1:length(year)){
  table <- bog_subset("TBH", lineages)
  table <- year_subset(year[i], table)
  abundance <- c()
  variance <- c()
  persistance <- c()
  
  for(j in 1:dim(table)[1]){
    row <- table[j, ]
    persistance[j] <- length(which(row > 0))/length(row)
    if(sum(row) > 0){
      variance[j] <- cv(as.numeric(row))
      abundance[j] <- sum(row[which(row > 0)])/length(which(row > 0))
    }else{
      variance[j] <- 0
      abundance[j] <- 0
    }
    
  }
  
  to.plot <- data.frame(abundance, persistance, variance)
  to.plot$lineage <- rownames(reduce_names(table))
  
  z <- ggplot(data = to.plot, aes(x = variance, y = abundance, fill = persistance)) + geom_point() + theme_bw() + labs(title = paste(year[i]), x = "Coefficient of Variance", y = "Mean Abundance when Present") + geom_label_repel(aes(label = lineage)) + scale_fill_gradient(low = "white", high = "lightgreen") + theme(legend.position = "none")
  
  assign(paste("p", i, sep = ""), z)
}

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS7.pdf", width = 3.3125*2, height = 7)
multiplot(p1, p2, p3, p4, cols = 2)
dev.off()


### Figure S8 - More examples of UniFrac vs Time

# Select lake and layer, then calculate UniFrac and plot vs time difference
NSH <- prune_samples(sampledata$Bog == "NS" & sampledata$Layer == "H", alldata)
x <- UniFrac(NSH, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p1 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "North Sparkling Hypolimnion", x = "Time Difference", y = "UniFrac Distance")

NSE <- prune_samples(sampledata$Bog == "NS" & sampledata$Layer == "E", alldata)
x <- UniFrac(NSE, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p2 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "North Sparkling Epilimnion", x = "Time Difference", y = "UniFrac Distance")

MAE <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer == "E", alldata)
x <- UniFrac(MAE, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p3 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Mary Lake Epilimnion", x = "Time Difference", y = "UniFrac Distance")

pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS8.1.pdf", width = 3.3125*2, height = 7)
multiplot(p1, p2, p3, cols = 1)
dev.off()

SSH <- prune_samples(sampledata$Bog == "SS" & sampledata$Layer == "H", alldata)
x <- UniFrac(SSH, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p1 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "South Sparkling Hypolimnion", x = "Time Difference", y = "UniFrac Distance")

SSE <- prune_samples(sampledata$Bog == "SS" & sampledata$Layer == "E", alldata)
x <- UniFrac(SSE, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p2 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "South Sparkling Epilimnion", x = "Time Difference", y = "UniFrac Distance")


pdf(file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS8.2.pdf", width = 3.3125*2, height = 7/3*2)
multiplot(p1, p2, cols = 1)
dev.off()