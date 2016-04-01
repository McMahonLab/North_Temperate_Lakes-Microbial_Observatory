library(OTUtable)
library(phyloseq)
library(ggplot2)
library(ape)
library(vegan)
library(grid)

data(otu_table)
data(taxonomy)

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
