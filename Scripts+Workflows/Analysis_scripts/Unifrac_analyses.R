library(OTUtable)
library(phyloseq)
library(ggplot2)
library(ape)
library(vegan)
library(grid)
library(reshape2)
path2repo <- "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/"
data(otu_table)
data(taxonomy)
data(metadata)
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


# Unfortunately, the tree of all OTUs is too large for ape to write to file (recursion errors)
# A new tree has to be made each time the script is run. This takes approx. 15 min for our 6,208 sequences. I'd recommend making yourself a cup of tea during this time.
seqs <- read.dna("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_repseqs_07Jul15.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

# Make OTU table, taxonomy, and sampledata datasets
OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(taxonomy))

sampledata <- sample_data(data.frame(Bog = substr(colnames(otu_table), start = 1, stop = 2), Layer = substr(colnames(otu_table), start = 3, stop = 3), Year = substr(colnames(otu_table), start = 9, stop = 10), row.names = colnames(otu_table), stringsAsfactors = F))                                                                                                                                          

alldata <- phyloseq(OTU, TAX, sampledata, bogtree)
years <- c("05", "07", "08", "09")

colors <- c("#a6cee3", "#1f78b4", "#33a02c", "#b2df8a")
# Make UniFrac PCoA of TBH
TBH <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "H", alldata)
TBH_year <- factor(substr(sample_names(TBH), start = 9, stop = 10), levels = years)

top <- names(sort(taxa_sums(TBH), TRUE)[1:500])
TBH_abun = prune_taxa(top, TBH)

x <- UniFrac(TBH_abun, weighted = T, normalize = T)
pcoa <- betadisper(x, TBH_year)
scores <- scores(pcoa)
# Locate centroids
TBHcentroids <- scores$centroids
TBHcentroids <- as.data.frame(TBHcentroids)
TBHcentroids$Year <- factor(years, level = years)

plot.pcoa <- data.frame(scores$sites, TBH_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

#Save for plotting
p1 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + geom_point(data=TBHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Trout Bog") + coord_cartesian(xlim = c(-0.3, 0.3), ylim = c(-0.4, 0.4)) + scale_color_manual(values = colors)

#Plot distance from centroids boxplots

TBH_distance <- pcoa$distance

# Calculate PERMADISP - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
adonis(x ~ Year, as(sample_data(TBH), "data.frame"))
# r2 0.34447
# Pr 0.001

# Calculate distance between centroids

# 05 to 07
sqrt((TBHcentroids[1,1] - TBHcentroids[2,1]) + (TBHcentroids[1,2] - TBHcentroids[2,2]))
# 07 to 08
sqrt((TBHcentroids[2,1] - TBHcentroids[3,1]) + (TBHcentroids[2,2] - TBHcentroids[3,2]))
# 08 to 09
sqrt((TBHcentroids[4,1] - TBHcentroids[3,1]) + (TBHcentroids[4,2] - TBHcentroids[3,2]))

#Repeat with North Sparkling
NSH <- prune_samples(sampledata$Bog == "NS" & sampledata$Layer == "H", alldata)
NSH_year <- factor(substr(sample_names(NSH), start = 9, stop = 10), levels = years)

top <- names(sort(taxa_sums(NSH), TRUE)[1:500])
NSH_abun = prune_taxa(top, NSH)

x <- UniFrac(NSH_abun, weighted = T, normalize = T)
pcoa <- betadisper(x, NSH_year)
scores <- scores(pcoa)
# Locate centroids
NSHcentroids <- scores$centroids
NSHcentroids <- as.data.frame(NSHcentroids)
NSHcentroids$Year <- factor(c("07", "08", "09"), level = years)

plot.pcoa <- data.frame(scores$sites, NSH_year[which(is.na(NSH_year) == F)])
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

#Save for plotting
p2 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + geom_point(data=NSHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "North Sparkling Bog") + coord_cartesian(xlim = c(-0.3, 0.3), ylim = c(-0.4, 0.4)) + scale_color_manual(values = colors[2:4]) 

NSH_distance <- pcoa$distance

# Calculate PERMADISP - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
adonis(x ~ Year, as(sample_data(NSH), "data.frame"))
# r2 0.11
# Pr 0.001
# Calculate distance between centroids

# 07 to 08
sqrt((NSHcentroids[2,1] - NSHcentroids[1,1]) + (NSHcentroids[2,2] - NSHcentroids[1,2]))
# 08 to 09
sqrt((NSHcentroids[2,1] - NSHcentroids[3,1]) + (NSHcentroids[2,2] - NSHcentroids[3,2]))

#Repeat with South Sparkling
SSH <- prune_samples(sampledata$Bog == "SS" & sampledata$Layer == "H", alldata)
SSH_year <- factor(substr(sample_names(SSH), start = 9, stop = 10), levels = years)

top <- names(sort(taxa_sums(SSH), TRUE)[1:500])
SSH_abun = prune_taxa(top, SSH)

x <- UniFrac(SSH_abun, weighted = T, normalize = T)
pcoa <- betadisper(x, SSH_year)
scores <- scores(pcoa)
# Locate centroids
SSHcentroids <- scores$centroids
SSHcentroids <- as.data.frame(SSHcentroids)
SSHcentroids$Year <- factor(c("07", "08", "09"), level = years)

plot.pcoa <- data.frame(scores$sites, SSH_year[which(is.na(SSH_year) == F)])
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

#Save for plotting
p3 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + geom_point(data=SSHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "South Sparkling Bog") + coord_cartesian(xlim = c(-0.3, 0.3), ylim = c(-0.4, 0.4)) + scale_color_manual(values = colors[2:4])

SSH_distance <- pcoa$distance

# Calculate PERMADISP - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
adonis(x ~ Year, as(sample_data(SSH), "data.frame"))
# r2 0.20488
# Pr 0.001
# Calculate distance between centroids

# 07 to 08
sqrt((SSHcentroids[1,1] - SSHcentroids[2,1]) + (SSHcentroids[1,2] - SSHcentroids[2,2]))
# 08 to 09
sqrt((SSHcentroids[3,1] - SSHcentroids[2,1]) + (SSHcentroids[3,2] - SSHcentroids[2,2]))


# Make UniFrac PCoA of MAH
MAH <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer == "H", alldata)
MAH_year <- factor(substr(sample_names(MAH), start = 9, stop = 10), levels = years)

top <- names(sort(taxa_sums(MAH), TRUE)[1:500])
MAH_abun = prune_taxa(top, MAH)

x <- UniFrac(MAH_abun, weighted = T, normalize = T)
pcoa <- betadisper(x, MAH_year)
scores <- scores(pcoa)
# Locate centroids
MAHcentroids <- scores$centroids
MAHcentroids <- as.data.frame(MAHcentroids)
MAHcentroids$Year <- factor(years, level = years)

plot.pcoa <- data.frame(scores$sites, MAH_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

#Save for plotting
p4 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + geom_point(data=MAHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Mary Lake") + coord_cartesian(xlim = c(-0.3, 0.3), ylim = c(-0.4, 0.4)) + scale_color_manual(values = colors)

MAH_distance <- pcoa$distance

# Calculate PERMADISP - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
adonis(x ~ Year, as(sample_data(MAH), "data.frame"))
# r2 0.10167
# Pr 0.002

# Calculate distance between centroids

# 05 to 07
sqrt((MAHcentroids[1,1] - MAHcentroids[2,1]) + (MAHcentroids[1,2] - MAHcentroids[2,2]))
# 07 to 08
sqrt((MAHcentroids[3,1] - MAHcentroids[2,1]) + (MAHcentroids[3,2] - MAHcentroids[2,2]))
# 08 to 09
sqrt((MAHcentroids[4,1] - MAHcentroids[3,1]) + (MAHcentroids[4,2] - MAHcentroids[3,2]))

pdf(file = paste("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/unifrac_pcoa_by_year.pdf", sep = ""), width = 3 * 2, height = 6)
multiplot(p1, p2, p3, p4, cols=2)
dev.off()

distances <- c(TBH_distance, NSH_distance, SSH_distance, MAH_distance)
lakes <- c(rep("TB", length(TBH_distance)), rep("NS", length(NSH_distance)), rep("SS", length(SSH_distance)), rep("MA", length(MAH_distance)))
lakes <- factor(lakes, level = c("TB", "NS", "SS", "MA"))

dispersion <- data.frame(distances, lakes)

pdf(file = paste("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/unifrac_pcoa_dispersion.pdf", sep = ""), width = 5, height = 3)
ggplot(data = dispersion, aes(y = distances, x = lakes, fill = lakes)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 8), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + scale_x_discrete(labels=c("Trout","N. Sparkling","S. Sparkling","Mary")) + labs(x = NULL, y = "Dispersion")
dev.off()

# Create a representative sample for Mary Lake for each year and test its similarity in PCoA

MAHtable <- bog_subset("MAH", otu_table)
MAHtable05 <- year_subset("05", MAHtable)
MAHtable07 <- year_subset("07", MAHtable)
MAHtable08 <- year_subset("08", MAHtable)
MAHtable09 <- year_subset("09", MAHtable)

otu_table2 <- otu_table
otu_table2$MAH01REP05 <- rowSums(MAHtable05)/dim(MAHtable05)[2]
otu_table2$MAH01REP07 <- rowSums(MAHtable07)/dim(MAHtable07)[2]
otu_table2$MAH01REP08 <- rowSums(MAHtable08)/dim(MAHtable08)[2]
otu_table2$MAH01REP09 <- rowSums(MAHtable09)/dim(MAHtable09)[2]


OTU2 <- otu_table(as.matrix(otu_table2), taxa_are_rows = T)
sampledata2 <- sample_data(data.frame(Bog = substr(colnames(otu_table2), start = 1, stop = 2), Layer = substr(colnames(otu_table2), start = 3, stop = 3), Year = substr(colnames(otu_table2), start = 9, stop = 10), row.names = colnames(otu_table2), stringsAsfactors = F)) 
alldata_reps <- phyloseq(OTU2, TAX, sampledata2, bogtree)

MAH2 <- prune_samples(sampledata2$Bog == "MA" & sampledata2$Layer == "H", alldata_reps)
MAH2_year <- factor(substr(sample_names(MAH2), start = 9, stop = 10), levels = years)
MAH2_type <- substr(sample_names(MAH2), start = 6, stop = 8) == "REP"

top <- names(sort(taxa_sums(MAH2), TRUE)[1:500])
MAH2_abun = prune_taxa(top, MAH2)

x <- UniFrac(MAH2_abun, weighted = T, normalize = T)
pcoa <- betadisper(x, MAH2_year)
scores <- scores(pcoa)
# Locate centroids
MAH2centroids <- scores$centroids
MAH2centroids <- as.data.frame(MAH2centroids)
MAH2centroids$Year <- factor(years, level = years)

plot.pcoa <- data.frame(scores$sites, MAH2_year, as.factor(MAH2_type))
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year", "Sample")

#Save for plotting
ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year, shape = Sample, size = Sample)) + geom_point() + scale_shape_manual(values = c(16, 15)) + scale_size_manual(values = c(2, 3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + geom_point(data=MAH2centroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Mary Lake") + coord_cartesian(xlim = c(-0.15, 0.1), ylim = c(-0.15, 0.1)) + scale_color_manual(values = colors)

# Representative samples for each year look great! Now test these against other lakes over time.
TBH07 <- prune_samples(sampledata2$Year == "07" & sampledata2$Bog == "TB" & sampledata2$Layer == "H" | sampledata2$Year == "07" & substr(sample_names(sampledata2), start = 6, stop = 8) == "REP", alldata_reps)
TBH07_lake <- factor(substr(sample_names(TBH07), start = 1, stop = 2), levels = c("TB", "MA"))

top <- names(sort(taxa_sums(TBH07), TRUE)[1:500])
TBH07_abun = prune_taxa(top, TBH07)
x <- UniFrac(TBH07, weighted = T, normalize = T)
sim <- 1 - as.matrix(x)[1:length(TBH07_lake)-1, length(TBH07_lake)]
TBH07_date <- extract_date(names(sim))
# plot(TBH07_date[order(TBH07_date)], sim[order(TBH07_date)], type = "l")

plot.trend <- data.frame(TBH07_date, sim)
colnames(plot.trend) <- c("Date", "Distance")

TBHmat <- make_temp_matrix("TBH.....07", metadata)
TBHmat <- melt(TBHmat)
colnames(TBHmat) <- c("Depth", "Date", "Temperature")
TBHmat$Date <- as.Date(TBHmat$Date, format = "%Y-%m-%d")
TBHmat$Depth <- -TBHmat$Depth / 11 + 0.79

# Need to close polygons - add 0 or max values at top and bottom of graph
TBHmat <- TBHmat[which(is.na(TBHmat$Temperature) == F), ]
# For each date, add a hidden value outside of the plotting range
add <- TBHmat[which(TBHmat$Depth == 0.34), ]
add$Depth <- rep(-1, length(add$Depth))
add$Temperature <- rep(4, length(add$Temperature))
add2 <- add
add2$Depth <- rep(3, length(add2$Depth))
add2$Temperature <- rep(28, length(add2$Temperature))
TBHmat2 <- rbind(TBHmat, add, add2)


pdf(file = paste(path2repo, "TBH_v_MAH_unifrac.pdf", sep = ""), width = 4, height = 2.5)
ggplot() + stat_contour(data = TBHmat2, aes(y = Depth, x = Date, z = Temperature, fill = ..level..), geom = "polygon") + scale_fill_gradientn(colours = c("dodgerblue", "cyan", "green", "yellow", "red"), "Temp", limits = c(4, 28)) + geom_line(data = plot.trend, aes(x = Date, y = Distance), size = 1.5) + labs(y = "1 - UniFrac Distance", x = NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "dodgerblue3"), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title.x = element_text(size = 15, vjust = 0.2), axis.title.y = element_text(size = 12, vjust = 1.6), axis.text.y = element_text(colour = "black", size = 10)) + coord_cartesian(xlim = extract_date(c("TBH20Jun07", "TBH11Nov07")), ylim = c(0.3, 0.75))
dev.off()

#Make PCoA of TB vs MA
TBH_MAH <- prune_samples(sampledata2$Year == "07" & sampledata2$Bog == "TB" & sampledata2$Layer == "H" | sampledata2$Year == "07"  & sampledata2$Bog == "MA" & sampledata2$Layer == "H" &  substr(sample_names(sampledata2), start = 6, stop = 8) != "REP", alldata_reps)
TBH_MAH_lake <- factor(substr(sample_names(TBH_MAH), start = 1, stop = 2), levels = c("TB", "MA"))

top <- names(sort(taxa_sums(TBH_MAH), TRUE)[1:500])
TBH_MAHabun = prune_taxa(top, TBH_MAH)
x <- UniFrac(TBH_MAHabun , weighted = T, normalize = T)
pcoa <- betadisper(x, TBH_MAH_lake)

scores <- scores(pcoa)
TBH_MAH_dates <- extract_date(sample_names(TBH_MAH))

plot.pcoa <- data.frame(scores$sites, TBH_MAH_lake, TBH_MAH_dates, estimate_richness(TBH_MAH, measures = "Shannon")/log(length(which(taxa_sums(TBH_MAH) > 0))))
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Lake", "Date", "Evenness")
plot.pcoa <- plot.pcoa[order(plot.pcoa$Date),]

#Calculate weights of axes
pcoa$eig[1]/sum(pcoa$eig)
pcoa$eig[2]/sum(pcoa$eig)

pdf(file = paste(path2repo, "TBH_v_MAH_unifrac_pcoa.pdf", sep = ""), width = 4, height = 2.75)
ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Evenness, shape = Lake)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_color_gradientn(colors = c("firebrick3", "gold2", "springgreen4")) + geom_path(data = plot.pcoa[which(plot.pcoa$Lake == "TB"),], aes(x = PCoA1, y = PCoA2), size = 0.25) + geom_path(data = plot.pcoa[which(plot.pcoa$Lake == "MA"),], aes(x = PCoA1, y = PCoA2), size = 0.25) + labs(x = "PCoA1 (56%)", y = "PCoA2 (14%)")
dev.off()

# North Sparkling over time
NSH08 <- prune_samples(sampledata2$Year == "08" & sampledata2$Bog == "TB" & sampledata2$Layer == "H" | sampledata2$Year == "08" & substr(sample_names(sampledata2), start = 6, stop = 8) == "REP", alldata_reps)
NSH08_lake <- factor(substr(sample_names(NSH08), start = 1, stop = 2), levels = c("TB", "MA"))

top <- names(sort(taxa_sums(NSH08), TRUE)[1:500])
NSH08_abun = prune_taxa(top, NSH08)
x <- UniFrac(NSH08_abun, weighted = T, normalize = T)
sim <- 1 - as.matrix(x)[1:length(NSH08_lake)-1, length(NSH08_lake)]
NSH08_date <- extract_date(names(sim))
plot(NSH08_date[order(NSH08_date)], sim[order(NSH08_date)], type = "l")

plot.trend <- data.frame(NSH08_date, sim)
colnames(plot.trend) <- c("Date", "Distance")

NSHmat <- make_temp_matrix("NSH.....08", metadata)
# Remove columns with NAs
NSHmat <- NSHmat[,c(1:14, 16:30, 32:54)]
NSHmat <- melt(NSHmat)
colnames(NSHmat) <- c("Depth", "Date", "Temperature")
NSHmat$Date <- as.Date(NSHmat$Date, format = "%Y-%m-%d")
NSHmat$Depth <- -NSHmat$Depth / 40 + 0.66
# Need to close polygons - add 0 or max values at top and bottom of graph
NSHmat <- NSHmat[which(is.na(NSHmat$Temperature) == F), ]
# For each date, add a hidden -1 value
add <- NSHmat[which(NSHmat$Depth == 0.66),]
add$Depth <- rep(0, length(add$Depth))
add$Temperature <- rep(4, length(add$Temperature))
add2 <- add
add2$Depth <- rep(3, length(add2$Depth))
add2$Temperature <- rep(28, length(add2$Temperature))
NSHmat2 <- rbind(NSHmat, add, add2)
# Remove dates with missing depth measurements

#pdf(file = paste(path2repo, "NSH_v_MAH_unifrac.pdf", sep = ""), width = 3.3125, height = 2.5)
ggplot() + stat_contour(data = NSHmat2, aes(y = Depth, x = Date, z = Temperature, fill = ..level..), geom = "polygon") + coord_cartesian(xlim = extract_date(c("NSH25May08", "NSH05Oct08")), ylim = c(0.32, 0.65)) + scale_fill_gradientn(colours = c("dodgerblue", "cyan", "green", "yellow", "red"), "Temp", limits = c(4, 28)) + geom_line(data = plot.trend, aes(x = Date, y = Distance), size = 1.5) + labs(y = "UniFrac Distance", x = NULL) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill="dodgerblue3"), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title.x = element_text(size = 15, vjust = 0.2), axis.title.y = element_text(size = 12, vjust = 1.6), axis.text.y = element_text(colour = "black", size = 10)) 
#dev.off()
#This plot doesn't look as nice because our metadata and DNA dates don't line up. Replace with a PCoA of TB vs MA?



# Use a GLM to test trends in similarity between lakes

hypo <- prune_samples(sampledata2$Layer == "H", alldata_reps)
top <- names(sort(taxa_sums(hypo), TRUE)[1:500])
hypo_abun = prune_taxa(top, hypo)
x <- UniFrac(hypo_abun, weighted = T, normalize = T)
x <- as.matrix(x)

sim <- c()
for(i in 1:dim(x)[1]){
  sampleyear <- substr(rownames(x)[i], start = 9, stop = 10)
  if (sampleyear == "05"){
    sim[i] <- x[i, 691]
  }else if (sampleyear == "07"){
    sim[i] <- x[i, 692]
  }else if (sampleyear == "08"){
    sim[i] <- x[i, 693]
  }else if (sampleyear == "09"){
    sim[i] <- x[i, 694]
  }
}
sim <- sim[1:690]

samplelakes <- substr(rownames(x)[1:690], start = 1, stop = 2)
sampleyear <- substr(rownames(x)[1:690], start = 9, stop = 10)
sampledate <- extract_date(rownames(x)[1:690])

glm.trend <- data.frame(rownames(x)[1:690], sim, samplelakes, sampleyear, sampledate)
colnames(glm.trend) <- c("Sample", "UniFrac", "Lake", "Year", "Date")

# Make a variable of continous mixing


cont.mixes <- c()
for(i in 1:dim(glm.trend)[1]){
  find <- grep(substr(glm.trend$Sample[i], start = 1, stop = 10), metadata$Sample_Name)
  if(length(find) > 0){
    samples <- metadata[find,]
    top <- samples$DO[which(samples$Depth == "0.5")]
    if(length(top) == 0){
      top <- samples$DO[which(samples$Depth == "1")]
    }
    bottom <- min(samples$DO)
    cont.mixes[i] <- top - bottom
  }
}

glm.trend$Mixing.cont <- cont.mixes
glm.trend$Mixing.bin <- glm.trend$Mixing.cont < 4

#Julian Date variables

julian <- c()
for(i in 1:dim(glm.trend)[1]){
  if(glm.trend$Year[i] == "05"){
    julian[i] <- as.numeric(glm.trend$Date[i] - extract_date(c("TBH01Jun05")))
  }else if(glm.trend$Year[i] == "07"){
    julian[i] <- as.numeric(glm.trend$Date[i] - extract_date(c("TBH01Jun07")))
  }else if(glm.trend$Year[i] == "08"){
    julian[i] <- as.numeric(glm.trend$Date[i] - extract_date(c("TBH01Jun08")))
  }else if(glm.trend$Year[i] == "09"){
    julian[i] <- as.numeric(glm.trend$Date[i] - extract_date(c("TBH01Jun09")))
  }
}

glm.trend$DayNum <- julian

#Remove Mary Lake - this comparison does not make sense
glm.trend <- glm.trend[which(glm.trend$Lake != "MA"),]

#hist(glm.trend$UniFrac)
#Distances look normal-ish

# Run linear model

# Run model without interaction to see main effects
model <- lm(UniFrac ~ Year + Lake + Mixing.bin + DayNum, data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model)

library(lme4)

model1 <- lmer(UniFrac ~ Year  + Mixing.bin + DayNum + (1|Lake), data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model1)

model2 <- lmer(UniFrac ~ Year  + Mixing.bin + Lake:DayNum + (1|Lake), data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model2)

model3 <- lmer(UniFrac ~ Mixing.bin + Lake:DayNum + (1|Lake) + (1|Year), data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model3)

model4 <- lmer(UniFrac ~ Mixing.bin + Lake:DayNum + (1|Lake:Year), data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model4)

model5 <- lmer(UniFrac ~ Mixing.cont + Lake:DayNum + (1|Lake:Year), data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model5)

model6 <- lm(UniFrac ~ Mixing.bin + Lake:DayNum:Year , data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model6)

model7 <- lm(UniFrac ~  Lake:DayNum:Year , data=glm.trend[which(glm.trend$DayNum >=0 & glm.trend$Mixing.bin == FALSE),])
summary(model7)

model8 <- lm(UniFrac ~  Lake:Mixing.cont:Year , data=glm.trend[which(glm.trend$DayNum >=0 ),])
summary(model8)

plot(glm.trend$Mixing.cont[which(glm.trend$Lake == "TB")] ~ glm.trend$DayNum[which(glm.trend$Lake == "TB")], col = glm.trend$Year[which(glm.trend$Lake == "TB")])

model <- lm(UniFrac ~ Year + Lake + Mixing.cont + DayNum, data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model)

#binary mixing seems to be the stronger predictor
#Is there a different effect of day number when mixing is true vs false?
model <- lm(UniFrac ~ Year + Lake + DayNum:Mixing.bin, data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model)
#Yes

#Is there a different effect of day number when mixing is false in each lake?
model <- lm(UniFrac ~ Year + Lake:DayNum:Mixing.bin, data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model)
#Unclear. Reduce the number of interactions.

#Is there a different effect of day number each lake?
model1 <- lm(UniFrac ~ Year + Lake*DayNum + Mixing.bin, data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model1)
#Need help interpreting this one

#Is there a different effect of mixing in each lake?
model2 <- lm(UniFrac ~ Year + DayNum + Lake:Mixing.bin, data=glm.trend[which(glm.trend$DayNum >=0),])
summary(model2)
#HK, FB, and CB all have a different effect - the most polymictic and meromictic? Maybe WS was essentially dimictic the one year tested.

# Check the residuals
plot(resid(model1) ~ fitted(model1))
plot(resid(model2) ~ fitted(model2))

#Model 1 actually looks nicer

# The model tells me that day number effects connectivity differently depending on how strongly stratified the lake is
# Check the predicted values for how day number affects different lakes
plot(fitted(model1)[which(glm.trend$DayNum >= 0 & glm.trend$Lake == "TB" & glm.trend$Mixing.bin  == F )], glm.trend$UniFrac[which(glm.trend$DayNum >= 0 & glm.trend$Lake == "TB" & glm.trend$Mixing.bin  == F )])

plot(fitted(model2)[which(glm.trend$DayNum >= 0 & glm.trend$Lake == "TB" & glm.trend$Mixing.bin  == F )], glm.trend$UniFrac[which(glm.trend$DayNum >= 0 & glm.trend$Lake == "TB" & glm.trend$Mixing.bin  == F )])
# Association between real and predicted values is not great

plot(glm.trend$DayNum[which(glm.trend$Lake == "TB" & glm.trend$Mixing.bin  == F )], glm.trend$UniFrac[which(glm.trend$Lake == "TB" & glm.trend$Mixing.bin  == F )])
#Negative linear decrease indicates increasing similarity over time


plot(glm.trend$DayNum[which(glm.trend$Lake == "HK" & glm.trend$Mixing.bin  == F )], glm.trend$UniFrac[which(glm.trend$Lake == "HK" & glm.trend$Mixing.bin  == F )])
#No trend

plot(glm.trend$DayNum[which(glm.trend$Lake == "CB" & glm.trend$Mixing.bin  == F )], glm.trend$UniFrac[which(glm.trend$Lake == "CB" & glm.trend$Mixing.bin  == F )])
#Slight negative linear decrease indicates increasing similarity over time

plot(glm.trend$DayNum[which(glm.trend$Lake == "SS" & glm.trend$Mixing.bin  == F )], glm.trend$UniFrac[which(glm.trend$Lake == "SS" & glm.trend$Mixing.bin  == F )])
#Slight negative linear decrease to no trend

plot(glm.trend$DayNum[which(glm.trend$Lake == "NS" & glm.trend$Mixing.bin  == F )], glm.trend$UniFrac[which(glm.trend$Lake == "NS" & glm.trend$Mixing.bin  == F )])
# Unclear

plot(glm.trend$DayNum[which(glm.trend$Lake == "WS" & glm.trend$Mixing.bin  == F )], glm.trend$UniFrac[which(glm.trend$Lake == "WS" & glm.trend$Mixing.bin  == F )])
#No trend
# where do unifrac and bray curtis disagree?

bc <- distance(hypo_abun, "bray")
bc <- as.matrix(bc)

sim <- c()
for(i in 1:dim(bc)[1]){
  sampleyear <- substr(rownames(bc)[i], start = 9, stop = 10)
  if (sampleyear == "05"){
    sim[i] <- x[i, 691]
  }else if (sampleyear == "07"){
    sim[i] <- x[i, 692]
  }else if (sampleyear == "08"){
    sim[i] <- x[i, 693]
  }else if (sampleyear == "09"){
    sim[i] <- x[i, 694]
  }
}
sim <- sim[1:690]
sim <- sim[which(samplelakes != "MA")]

glm.trend$BrayCurtis <- sim 

#Try same model with BC
model7.2 <- lm(BrayCurtis ~  Lake:DayNum:Year , data=glm.trend[which(glm.trend$DayNum >=0 & glm.trend$Mixing.bin == FALSE),])
summary(model7.2)

model8.2 <- lm(BrayCurtis ~  Lake:Mixing.cont:Year , data=glm.trend[which(glm.trend$DayNum >=0 ),])
summary(model8.2)

#Same results - is it the abundance threshold?

bc <- distance(hypo, "bray")
bc <- as.matrix(bc)

sim <- c()
for(i in 1:dim(bc)[1]){
  sampleyear <- substr(rownames(bc)[i], start = 9, stop = 10)
  if (sampleyear == "05"){
    sim[i] <- x[i, 691]
  }else if (sampleyear == "07"){
    sim[i] <- x[i, 692]
  }else if (sampleyear == "08"){
    sim[i] <- x[i, 693]
  }else if (sampleyear == "09"){
    sim[i] <- x[i, 694]
  }
}
sim <- sim[1:690]
sim <- sim[which(samplelakes != "MA")]

glm.trend$BrayCurtis <- sim 

simple <- lm(BrayCurtis ~ Lake:DayNum:Year, data = glm.trend[which(glm.trend$Mixing.bin == F),])

#Try re-running UniFrac comparing to all Mary Lake samples instead of the representative

MAsamples <- x[which(substr(rownames(x), start = 1, stop = 2) != "MA"),which(substr(rownames(x), start = 1, stop = 2) == "MA" & substr(rownames(x), start = 6, stop = 8) != "REP")]
mean.Unifrac <- rowSums(MAsamples)/dim(MAsamples)[2]

glm.trend$meanUnifrac <- mean.Unifrac

model7.3 <- lm(meanUnifrac ~  Lake:DayNum:Year , data=glm.trend[which(glm.trend$DayNum >=0 & glm.trend$Mixing.bin == FALSE),])
summary(model7.3)

model8.3 <- lm(meanUnifrac ~  Lake:Mixing.cont:Year , data=glm.trend[which(glm.trend$DayNum >=0 ),])
summary(model8.3)

#Results still odd - add in the mean Bray Curtis and see if I get the same result as previously
MAsamples <- bc[which(substr(rownames(x), start = 1, stop = 2) != "MA"),which(substr(rownames(x), start = 1, stop = 2) == "MA" & substr(rownames(x), start = 6, stop = 8) != "REP")]
mean.bc <- rowSums(MAsamples)/dim(MAsamples)[2]

glm.trend$meanBrayCurtis <- mean.bc

model7.4 <- lm(meanBrayCurtis ~  Lake:DayNum:Year , data=glm.trend[which(glm.trend$DayNum >=0 & glm.trend$Mixing.bin == FALSE),])
summary(model7.4)

model8.4 <- lm(meanBrayCurtis ~  Lake:Mixing.cont:Year , data=glm.trend[which(glm.trend$DayNum >=0 ),])
summary(model8.4)

#Results are still off. Try removing the abundance thresholds.

bc <- distance(hypo, "bray")
bc <- as.matrix(bc)
x <- UniFrac(hypo, weighted = T, normalize = T)
x <- as.matrix(x)

#Start trying to find where bray curtis and unifrac disagree

plot(glm.trend$meanBrayCurtis, glm.trend$meanUnifrac)
summary(lm(glm.trend$meanBrayCurtis ~ glm.trend$meanUnifrac))

#Highly correlated - this means that the difference is the clade vs otu designation, not UniFrac vs Bray Curtis

clade_table <- combine_otus("Clade", otu_table, taxonomy)
hypo_clade_table <- bog_subset("..H", clade_table)
bc.clade <- vegdist(t(hypo_clade_table), method = "bray")
bc.clade <- as.matrix(bc.clade)
bc.clade <- bc.clade[which(substr(rownames(bc.clade), start = 1, stop = 2) != "MA"), which(substr(rownames(bc.clade), start = 1, stop = 2) == "MA" )]

glm.trend$cladeBrayCurtis <- rowSums(bc.clade)/dim(bc.clade)[2]

plot(glm.trend$cladeBrayCurtis, glm.trend$meanUnifrac)
summary(lm(glm.trend$meanBrayCurtis ~ glm.trend$meanUnifrac))

#Also good correlation
plot(glm.trend$cladeBrayCurtis, glm.trend$meanBrayCurtis, col = glm.trend$Lake)

model7.4 <- lm(cladeBrayCurtis ~  Lake:DayNum:Year , data=glm.trend[which(glm.trend$DayNum >=0 & glm.trend$Mixing.bin == FALSE),])
summary(model7.4)
