#manuscript v0 figures

###########
# set up
library(OTUtable)
library(reshape2)
library(ggplot2)
library(ape)
library(phyloseq)
library(vegan)
library(VennDiagram)
library(exactRankTests)
library(cowplot)
library(raster)
library(ggrepel)
library(scales)
#library(indicspecies)

data(metadata)
data(otu_table)
data(taxonomy)
seqs <- read.dna("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_repseqs_07Jul15.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

# #find average ph by lake in 2007 - for table 1. Input file is 2007 metadata file from my D drive, organized bog data/
# datasheet <- read.csv(file.choose(), header = T)
# phdata <- datasheet[which(is.na(datasheet$pH) == F), ]
# phdata$LakeLayer <- substr(phdata$Sample_Name, start = 1, stop = 3)
# red_phdata <- phdata[, c(8, 11)]
# attach(red_phdata)
# meandata <- aggregate(red_phdata, by = list(pH, LakeLayer), FUN = mean)
# 
# library(plyr)
# ddply(red_phdata, .(LakeLayer), summarize, mean_value = mean(pH))


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

# For a supplmentary table, report the mean and standard deviation of each lake-layer
mean(epi.data$epi.obs[which(epi.data$epi.lakes == "CB")]) #129
sd(epi.data$epi.obs[which(epi.data$epi.lakes == "CB")]) #28
mean(epi.data$epi.obs[which(epi.data$epi.lakes == "FB")]) #109
sd(epi.data$epi.obs[which(epi.data$epi.lakes == "FB")]) #32
mean(epi.data$epi.obs[which(epi.data$epi.lakes == "WS")]) #150
sd(epi.data$epi.obs[which(epi.data$epi.lakes == "WS")]) #45
mean(epi.data$epi.obs[which(epi.data$epi.lakes == "NS")]) #143
sd(epi.data$epi.obs[which(epi.data$epi.lakes == "NS")]) #33
mean(epi.data$epi.obs[which(epi.data$epi.lakes == "TB")]) #148
sd(epi.data$epi.obs[which(epi.data$epi.lakes == "TB")]) #38
mean(epi.data$epi.obs[which(epi.data$epi.lakes == "SS")]) #191
sd(epi.data$epi.obs[which(epi.data$epi.lakes == "SS")]) #57
mean(epi.data$epi.obs[which(epi.data$epi.lakes == "HK")]) #199
sd(epi.data$epi.obs[which(epi.data$epi.lakes == "HK")]) #67
mean(epi.data$epi.obs[which(epi.data$epi.lakes == "MA")]) #259
sd(epi.data$epi.obs[which(epi.data$epi.lakes == "MA")]) #67


mean(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "CB")]) #148
sd(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "CB")]) #31
mean(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "FB")]) #145
sd(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "FB")]) #57
mean(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "WS")]) #182
sd(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "WS")]) #56
mean(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "NS")]) #178
sd(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "NS")]) #40
mean(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "TB")]) #186
sd(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "TB")]) #38
mean(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "SS")]) #191
sd(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "SS")]) #54
mean(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "HK")]) #397
sd(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "HK")]) #124
mean(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "MA")]) #477
sd(hypo.data$hypo.obs[which(hypo.data$hypo.lakes == "MA")]) # 110

# plot 1A
fig1a <- ggplot(data = epi.data, aes(y = epi.obs, x = epi.lakes, fill = epi.lakes)) + geom_boxplot() + labs(y = "Observed Richness", x = NULL) + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_blank()) + labs(title = "Epilimnion") + ylim(0, 650) + theme(legend.position = "none") + background_grid(major = "xy", minor = "none")

# plot 1B
fig1b <- ggplot(data = hypo.data, aes(y = hypo.obs, x = hypo.lakes, fill = hypo.lakes)) + geom_boxplot() + labs(y="Observed Richness", x = NULL) + scale_fill_brewer(palette = "Paired") + labs(title = "Hypolimnion") + ylim(0, 650) + theme(legend.position = "none") + background_grid(major = "xy", minor = "none")

fig1 <- plot_grid(fig1a, fig1b, nrow = 2, labels = c("A", "B"), align = "w")
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/Fig1.pdf", fig1, base_aspect_ratio = 1.3, base_height = 4.75)


# Check significance using the Wilcoxon Rank Test on medians (not output as pdf, indicated as symbols in Illustrator)
pairwise.wilcox.test(epi.data$epi.obs, epi.data$epi.lakes, p.adjust.method = "bonferroni")

pairwise.wilcox.test(hypo.data$hypo.obs, hypo.data$hypo.lakes, p.adjust.method = "bonferroni")

##########
#Figure 2

# # Set up phyloseq object to run UniFrac on
OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(taxonomy))
sampledata <- sample_data(data.frame(Bog = substr(colnames(otu_table), start = 1, stop = 2), Layer = substr(colnames(otu_table), start = 3, stop = 3), Year = substr(colnames(otu_table), start = 9, stop = 10), row.names = colnames(otu_table), stringsAsfactors = F)) 
alldata <- phyloseq(OTU, TAX, sampledata, bogtree)

# Separate epilimnion and hypolimnion samples
epi <- prune_samples(sampledata$Layer == "E", alldata)
hypo <- prune_samples(sampledata$Layer == "H", alldata)

# Analyze and plot epilimnion points
#epi <- prune_taxa(taxa_sums(alldata) > 1000, alldata)
x <- UniFrac(epi, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 1, stop = 3))
scores <- scores(pcoa)
lakes <- factor(substr(rownames(scores$sites), start = 1, stop = 2), levels = c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA"))
regime <- c()
regime[which(lakes == "CB" | lakes == "FB" | lakes == "WS")] <- "polymictic"
regime[which(lakes == "TB" | lakes == "NS" | lakes == "SS")] <- "dimictic"
regime[which(lakes == "HK" | lakes == "MA")] <- "meromictic"

plot.pcoa <- data.frame(scores$sites, lakes)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Lake")

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(plot.pcoa$Lake)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(plot.pcoa[plot.pcoa$Lake == g, ],
                                                   veganCovEllipse(cov.wt(cbind(PCoA1,PCoA2), wt = rep(1/length(PCoA1), length(PCoA1)))$cov , center=c(mean(PCoA1), mean(PCoA2)))))
                                , group = g))
}

colnames(df_ell) <- c("PCoA1", "PCoA2", "Lake")

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

fig2a <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1*100, "% variability)", sep = ""), y = paste("PCoA2 (", axis2*100, "% variability)")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + geom_path(data = df_ell, aes(x = PCoA1, y = PCoA2, colour = Lake), size = 1, linetype = 1) + theme(legend.title = element_blank())

adonis(x ~ lakes, as(sample_data(epi), "data.frame"))
# r^2 = 0.34498, p = 0.001
# Calculate PERMANOVA by mixing regime
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
plot.pcoa <- data.frame(scores$sites, lakes)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2",  "Lake")
# Calculate percent variation explained of each axis
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

df_ell <- data.frame()
for(g in levels(plot.pcoa$Lake)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(plot.pcoa[plot.pcoa$Lake == g, ],
                                                   veganCovEllipse(cov.wt(cbind(PCoA1,PCoA2), wt = rep(1/length(PCoA1), length(PCoA1)))$cov , center=c(mean(PCoA1), mean(PCoA2)))))
                                , group = g))
}
colnames(df_ell) <- c("PCoA1", "PCoA2", "Lake")

fig2b <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1*100, "% variability)", sep = ""), y = paste("PCoA2 (", axis2*100, "% variability)")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + geom_path(data = df_ell, aes(x = PCoA1, y = PCoA2, colour = Lake), size = 1, linetype = 1) + theme(legend.position = "none")


# Calculate significance level of clustering by group
adonis(x ~ lakes, as(sample_data(hypo), "data.frame"))
# r^2 = 0.48651, p = 0.001
# Calculate PERMANOVA by mixing regime
adonis(x ~ regime, as(cbind(sample_data(hypo), regime), "data.frame"))
# r^2 = 0.21907, p = 0.001

legend <- get_legend(fig2a)
fig2a <- fig2a + theme(legend.position = "none")

fig2 <- plot_grid(fig2a, fig2b, nrow = 2, labels = c("A", "B"), align = "w")
fig2 <- plot_grid(fig2, legend, nrow = 1, rel_widths = c(1, 0.1))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/Fig2.pdf", fig2, base_aspect_ratio = 1.2, base_height = 6)

##########
#Figure 3
# 3 PCoAs by lake by year - TB, SS, and MA
# Then boxplot of dispersion by lake and layer

years <- c("05", "07", "08", "09")

colors <- c("#a6cee3", "#1f78b4", "#33a02c", "#b2df8a")
# Make UniFrac PCoA of TBH
TBH <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "H", alldata)
TBH_year <- factor(substr(sample_names(TBH), start = 9, stop = 10), levels = years)

x <- UniFrac(TBH, weighted = T, normalize = T)
pcoa <- betadisper(x, TBH_year)
scores <- scores(pcoa)
# Locate centroids
TBHcentroids <- scores$centroids
TBHcentroids <- as.data.frame(TBHcentroids)
TBHcentroids$Year <- factor(years, level = years)

plot.pcoa <- data.frame(scores$sites, TBH_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

fig3a <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + geom_point(data=TBHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Trout Bog", x = paste("PCoA1 (", axis1*100, "%)", sep = ""), y = paste("PCoA2 (", axis2*100, "%)"), subtitle = "r2 = 0.36, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors) + theme(plot.subtitle = element_text(hjust = 0.5))
fig3_legend <- get_legend(fig3a)
fig3a <- fig3a + theme(legend.position = "none")

# Calculate PERMANOVA - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
adonis(x ~ Year, as(sample_data(TBH), "data.frame"))
# r2 0.35531
# Pr 0.001

# PCoA of South Sparkling
SSH <- prune_samples(sampledata$Bog == "SS" & sampledata$Layer == "H", alldata)
SSH_year <- factor(substr(sample_names(SSH), start = 9, stop = 10), levels = years)

x <- UniFrac(SSH, weighted = T, normalize = T)
pcoa <- betadisper(x, SSH_year)
scores <- scores(pcoa)
# Locate centroids
SSHcentroids <- scores$centroids
SSHcentroids <- as.data.frame(SSHcentroids)
SSHcentroids$Year <- factor(years[2:4], level = years)

plot.pcoa <- data.frame(scores$sites, SSH_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

fig3b <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + geom_point(data=SSHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "South Sparkling Bog", x = paste("PCoA1 (", axis1*100, "%)", sep = ""), y = paste("PCoA2 (", axis2*100, " %)"), subtitle = "r2 = 0.20, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors[2:4])  + theme(legend.position = "none", plot.subtitle = element_text(hjust = 0.5))

# Calculate PERMANOVA - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
adonis(x ~ Year, as(sample_data(SSH), "data.frame"))
# r2 0.20469
# Pr 0.001

# PCoA of Mary Lake
MAH <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer == "H", alldata)
MAH_year <- factor(substr(sample_names(MAH), start = 9, stop = 10), levels = years)

x <- UniFrac(MAH, weighted = T, normalize = T)
pcoa <- betadisper(x, MAH_year)
scores <- scores(pcoa)
# Locate centroids
MAHcentroids <- scores$centroids
MAHcentroids <- as.data.frame(MAHcentroids)
MAHcentroids$Year <- factor(years, level = years)

plot.pcoa <- data.frame(scores$sites, MAH_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

fig3c <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=MAHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Mary Lake", x = paste("PCoA1 (", axis1*100, "%)", sep = ""), y = paste("PCoA2 (", axis2*100, "%)"), subtitle = "r2 = 0.09, p = 0.002") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)

# Calculate PERMANOVA - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
adonis(x ~ Year, as(sample_data(MAH), "data.frame"))
# r2 0.09169
# Pr 0.002

#Panel 2

# Use weighted UniFrac matrix to measure mean pairwise distance within sites and years
# Make scatterplots with mean as line - geom_dotplot

beta_diversity <- UniFrac(alldata, weighted = T, normalize = T)
site <- c("CBE", "FBE", "WSE", "NSE", "TBE", "SSE", "HKE", "MAE", "CBH", "FBH", "WSH", "NSH", "TBH", "SSH", "HKH", "MAH")
melted_diversity <- melt(as.matrix(beta_diversity))
melted_diversity$Site1 <- substr(melted_diversity$Var1, start = 1, stop = 3)
melted_diversity$Site2 <- substr(melted_diversity$Var2, start = 1, stop = 3)
within_site <- melted_diversity[which(melted_diversity$Site1 == melted_diversity$Site2), ]

#Remove anything with the letter U
within_site <- within_site[grep("U", within_site$Site1, invert = T), ]
within_site$Site1 <- factor(within_site$Site1, levels = rev(c("CBE","CBH", "FBE", "FBH", "WSE", "WSH", "NSE", "NSH", "TBE", "TBH", "SSE", "SSH", "HKE", "HKH", "MAE", "MAH")))
pal <- rev(c("#a6cee3", "#a6cee3", "#1f78b4", "#1f78b4", "#b2df8a", "#b2df8a", "#33a02c", "#33a02c", "#fb9a99", "#fb9a99", "#e31a1c", "#e31a1c", "#fdbf6f", "#fdbf6f", "#ff7f00", "#ff7f00"))
means <- c()
for(i in 1:length(site)){
  values <- within_site$value[which(within_site$Site1 == site[i])]
  means[i] <- mean(values)
}
within_site2 <- data.frame(site, means)
colnames(within_site2) <- c("Site1", "value")

fig3d <- ggplot(data = within_site, aes(x = Site1, y = value, fill = Site1)) + geom_violin() + coord_flip() + theme(legend.position = "none") + scale_fill_manual(values = pal) + labs(y = "Pairwise Weighted UniFrac Distance", x = NULL) + geom_point(data = within_site2, aes(y = value, x = Site1), pch = 16, size = 3)

#Test pairwise significance
pairwise.wilcox.test(within_site$value, within_site$Site1, p.adjust.method = "bonferroni")

fig3_panel1 <- plot_grid(fig3a, fig3b, fig3c, nrow = 3, labels = c("A", "B", "C"))
fig3_panel1 <- plot_grid(fig3_panel1, fig3_legend, nrow = 1, rel_widths = c(1, 0.1))
fig3 <- plot_grid(fig3_panel1, fig3d, nrow = 1, rel_widths = c(1, 1), labels = c("", "D"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/Fig3.pdf", fig3, base_aspect_ratio = 1.5, base_height = 6)

###############
#Figure 4
#Venn Diagrams

#Part 1 - core analysis, in text

# Is there a natural breaking point?
FWcore <- rownames(otu_table)[which(rowSums(otu_table > 0) == 1387)]
FWcore <- taxonomy[match(FWcore, rownames(taxonomy)), ]

persistant_levels <- c()
for(i in 1:100){
  persistant_levels[i] <- length(which(rowSums(otu_table > 0) >= 1387*(i/100)))
}
barplot(persistant_levels)
# No natural breaking point

#Otu0097 is PnecC
#What about most samples?
FWcore95 <- rownames(otu_table)[which(rowSums(otu_table > 0) >= 1387*0.95)]
FWcore95 <- taxonomy[match(FWcore95, rownames(taxonomy)), ]
#76, 97, 813 == bacI-A1, PnecC, acI-B2
FWcore90 <- rownames(otu_table)[which(rowSums(otu_table > 0) >= 1387*0.90)]
FWcore90 <- taxonomy[match(FWcore90, rownames(taxonomy)), ]
# add 678 == LD28

epi <- bog_subset("..E", otu_table)
Epicore90 <- rownames(epi)[which(rowSums(epi > 0) >= 671*0.90)]
Epicore90 <- taxonomy[match(Epicore90, rownames(taxonomy)), ]
# 4, 76, 97, 184, 472, 522, 678, 813 == betI, bacI-A1, PnecC, acI-B3, Lhab-A4, alfI-A1, LD28, acI-B2
hypo <- bog_subset("..H", otu_table)
hypocore90 <- rownames(hypo)[which(rowSums(hypo > 0) >= 690*0.90)]
hypocore90 <- taxonomy[match(hypocore90, rownames(taxonomy)), ]
# 42, 53, 76, 97, 189, 813 == Rhodo, unclassified Verruco, bacI-A1, PnecC, acI-B2, acI-B2

#Repeat with tribes
tribe_table <- combine_otus("Tribe", otu_table, taxonomy)
FWcore <- rownames(tribe_table)[which(rowSums(tribe_table > 0) == 1387)]
FWcore <- taxonomy[match(FWcore, rownames(taxonomy)), ]
FWcore95 <- rownames(tribe_table)[which(rowSums(tribe_table > 0) >= 1387*0.95)]
FWcore95 <- taxonomy[match(FWcore95, rownames(taxonomy)), ]
FWcore90 <- rownames(tribe_table)[which(rowSums(tribe_table > 0) >= 1387*0.90)]
FWcore90 <- taxonomy[match(FWcore90, rownames(taxonomy)), ]

epi <- bog_subset("..E", tribe_table)
Epicore90 <- rownames(epi)[which(rowSums(epi > 0) >= 671*0.90)]
Epicore90 <- taxonomy[match(Epicore90, rownames(taxonomy)), ]
hypo <- bog_subset("..H", tribe_table)
hypocore90 <- rownames(hypo)[which(rowSums(hypo > 0) >= 690*0.90)]
hypocore90 <- taxonomy[match(hypocore90, rownames(taxonomy)), ]

# Part 2 - overlap by mixing regime

FBE.vector <- rowSums(bog_subset("FBE", otu_table))
CBE.vector <- rowSums(bog_subset("CBE", otu_table))
WSE.vector <- rowSums(bog_subset("WSE", otu_table))
NSE.vector <- rowSums(bog_subset("NSE", otu_table))
TBE.vector <- rowSums(bog_subset("TBE", otu_table))
SSE.vector <- rowSums(bog_subset("SSE", otu_table))
HKE.vector <- rowSums(bog_subset("HKE", otu_table))
MAE.vector <- rowSums(bog_subset("MAE", otu_table))

FBH.vector <- rowSums(bog_subset("FBH", otu_table))
CBH.vector <- rowSums(bog_subset("CBH", otu_table))
WSH.vector <- rowSums(bog_subset("WSH", otu_table))
NSH.vector <- rowSums(bog_subset("NSH", otu_table))
TBH.vector <- rowSums(bog_subset("TBH", otu_table))
SSH.vector <- rowSums(bog_subset("SSH", otu_table))
HKH.vector <- rowSums(bog_subset("HKH", otu_table))
MAH.vector <- rowSums(bog_subset("MAH", otu_table))


poly.epi <- FBE.vector > 0 & CBE.vector > 0 & WSE.vector > 0
poly.hypo <- FBH.vector > 0 & CBH.vector > 0 & WSH.vector > 0
di.epi <- NSE.vector > 0 & TBE.vector > 0 & SSE.vector > 0
di.hypo <- NSH.vector > 0 & TBH.vector > 0 & SSH.vector > 0
mero.epi <- HKE.vector > 0 & MAE.vector > 0
mero.hypo <- HKH.vector > 0 & MAH.vector > 0

grid.newpage()
v1 <- draw.triple.venn(area1 = length(which(poly.epi == T)), area2 = length(which(di.epi == T)), area3 = length(which(mero.epi == T)), n12 = length(which(poly.epi == T & di.epi == T)), n13 = length(which(poly.epi == T & mero.epi == T)), n23 = length(which(mero.epi == T & di.epi == T)), n123 = length(which(poly.epi == T & di.epi == T & mero.epi == T)), fill = c("#a6cee3", "#33a02c", "#fdbf6f"), alpha = 0.4)

grid.newpage()
v2 <- draw.triple.venn(area1 = length(which(poly.hypo == T)), area2 = length(which(di.hypo == T)), area3 = length(which(mero.hypo == T)), n12 = length(which(poly.hypo == T & di.hypo == T)), n13 = length(which(poly.hypo == T & mero.hypo == T)), n23 = length(which(mero.hypo == T & di.hypo == T)), n123 = length(which(poly.hypo == T & di.hypo == T & mero.hypo == T)), fill = c("#a6cee3", "#33a02c", "#fdbf6f"), alpha = 0.4)

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/Fig4a.pdf", width = 2, height = 2)
grid.draw(v1)
dev.off()

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/Fig4b.pdf", width = 2, height = 2)
grid.draw(v2)
dev.off()

##########
#### Figure 5
# Measure traits of lineages to show they are consistent
taxonomy2 <- clean_TaxAss_taxonomy("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bogs_reclassified_11Mar16.csv", otu_table, remove_bootstrap = T)

lineage_table <- combine_otus("Lineage", otu_table, taxonomy2)

interest <- c(";betII$", ";acI$", ";betI$", ";bacI$", ";gamI$", ";verI-A$", ";betIV$", ";acV$", ";gamIII$", ";alfI$", ";betIII$")

lineages <- lineage_table[grep(paste(interest, collapse = "|"), rownames(lineage_table)), ]
# I'm only showing a few lakes, but these can be adjusted easilty
lakes <- c("CBE", "TBE", "SSE", "MAE")
lakenames <- c("Crystal Bog", "Trout Bog", "South Sparkling Bog", "Mary Lake")

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
  
  z <- ggplot(data = to.plot, aes(x = variance, y = abundance/2500, fill = persistence)) + geom_point() + theme_bw() + labs(title = paste(lakenames[i]), x = "Variability (CV)", y = "Mean Abundance when Present") + geom_label_repel(aes(label = lineage)) + scale_fill_gradient(low = "white", high = "lightgreen") + theme(legend.title = element_blank())
  assign(paste("p", i, sep = ""), z)
}
fig5_legend <- get_legend(p2)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

fig5 <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = c("A", "B", "C", "D"))
fig5 <- plot_grid(fig5, fig5_legend, nrow = 1, rel_widths = c(1, 0.1))

save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/Fig5.pdf", fig5, base_aspect_ratio = 1.5, base_height = 6)

##########
#Supplemental figures
##########
#Figure S1
#Sampling frequency

site <- substr(colnames(otu_table), start = 1, stop = 3)
date <- extract_date(colnames(otu_table))
sampling_freq <- data.frame(colnames(otu_table), site, date)
colnames(sampling_freq) <- c("Sample", "Site", "Date")

figs1 <- ggplot(sampling_freq, aes(x = Date, y = Site)) + geom_point() + scale_x_date(date_breaks = "months", date_labels = "%b-%y") + theme(axis.text.x = element_text(angle = -90, hjust = 1))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS1.pdf", figs1, base_aspect_ratio = 2, base_height = 4)


#Figure S2
#Phylum rank abundance

phylum_table <- combine_otus("Phylum", otu_table, taxonomy2)
archaea <- grep("Archaea", rownames(phylum_table))

phylum_table <- reduce_names(phylum_table)

# Sum total observations of each phylum in the full dataset and sort by abundance
totals <- rowSums(phylum_table)
totals <- sort(totals)

# Remove the "p__" phylum designation from phylum names. This looks better for plotting.
get.names <- strsplit(names(totals), "p__")
phyla.names <- c()
for(i in 1:length(get.names)){
  phyla.names[i] <- get.names[[i]][2]
}
phyla.names[23] <- "unclassified Archaea"
phyla.names[30] <- "unclassified"
phyla.names[45] <- "unclassified Bacteria"

# Set up a dataframe for plotting in ggplot2. Set the phyla.names to factors, with levels in order of abundance.
phylum_totals <- data.frame(phyla.names, totals)
phylum_totals$phyla.names <- factor(phylum_totals$phyla.names, levels = phylum_totals$phyla.names[order(phylum_totals$totals)])

# Plot as a bar graph

figs2 <- ggplot(data=phylum_totals, aes(x = phyla.names, y = totals)) + geom_bar(stat = "identity", fill = "grey", colour="black") + labs(x = NULL, y = "Log of Observed Reads") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, colour = "black"), axis.title = element_text(size = 15, vjust=1.2), plot.title = element_text(size = 20), axis.text.y = element_text(colour="black")) + scale_y_log10(expand = c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS2.pdf", figs2, base_aspect_ratio = 2, base_height = 5)

###########
#Figure S3
# Mega PCoA!
x <- UniFrac(prune_samples(sampledata$Layer != "U", alldata), weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 1, stop = 3))
scores <- scores(pcoa)
layer <- factor(substr(rownames(scores$sites), start = 3, stop = 3), levels = c("E", "H"))
lakes <- factor(substr(rownames(scores$sites), start = 1, stop = 2), levels = c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA"))
regime <- c()
regime[which(lakes == "CB" | lakes == "FB" | lakes == "WS")] <- "polymictic"
regime[which(lakes == "TB" | lakes == "NS" | lakes == "SS")] <- "dimictic"
regime[which(lakes == "HK" | lakes == "MA")] <- "meromictic"

plot.pcoa <- data.frame(scores$sites, layer, lakes)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Layer", "Lake")

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

figS3 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake, shape = Layer)) + geom_point(size = 3, alpha = 1/2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), panel.border = element_rect(colour = "black", fill=NA, size=1)) + labs(title = "All Data", x = paste("PCoA1 (", axis1*100, "% variation explained)", sep = ""), y = paste("PCoA2 (", axis2*100, "% variation explained)")) + scale_shape_manual(values=c(21, 22, 23, 24))  + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired")
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS3.pdf", figS3, base_aspect_ratio = 2.5, base_height = 4)

##########

#Figure S4 - richness over time

# Trout Bog, 2007

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

figS4a <- ggplot() + geom_line(data = TB_richness, aes(x = hypo.date, y = hypo.rich), size = 1) + labs(title = "Trout Bog", x = NULL, y = "Observed Richness")

# Repeat with North Sparkling, 2008

hypo <- bog_subset(paste("NSH", sep = ""), otu_table)
hypo <- year_subset("08", hypo)
hypo.rich <- apply(hypo, 2, obs_richness)
hypo.date <- extract_date(colnames(hypo))

NS_richness <- data.frame(hypo.date, hypo.rich)


figS4b <- ggplot() + geom_line(data = NS_richness, aes(x = hypo.date, y = hypo.rich), size = 1) + labs(title = "Trout Bog", x = NULL, y = "Observed Richness")


figS4 <- plot_grid(figS4a, figS4b, align = "h", nrow = 2, labels = c("A", "B"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS4.pdf", figS4, base_aspect_ratio = 1.5, base_height = 5)
plot_column(make_temp_matrix("TBE.....07", metadata), "Trout Bog 2007")
plot_column(make_temp_matrix("NSE.....08", metadata), "North Sparkling Bog 2008")

#add the water column plots in illustrator


##########
#Figure S5 - Run the UniFrac by lake by layer
CB <- prune_samples(sampledata$Bog == "CB" & sampledata$Layer != "U", alldata)
x <- UniFrac(CB, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
scores <- scores(pcoa)
layer <- substr(labels(x), start = 3, stop = 3)
plot.pcoa <- data.frame(scores$sites, layer)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Layer")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
plot1 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Layer)) + geom_point(size=2, alpha = 0.5) + labs(title = "Crystal Bog, r2 = 0.03", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_manual(values = c("#d8b365", "#5ab4ac"))
#plot_legend <- get_legend(plot1)
plot1 <- plot1 + theme(legend.position = "none")
#adonis(x ~ layer, as(sample_data(CB), "data.frame")) #0.03235 p = 0.003

FB <- prune_samples(sampledata$Bog == "FB" & sampledata$Layer != "U", alldata)
x <- UniFrac(FB, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
scores <- scores(pcoa)
layer <- substr(labels(x), start = 3, stop = 3)
plot.pcoa <- data.frame(scores$sites, layer)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Layer")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
plot2 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Layer)) + geom_point(size=2, alpha = 0.5) + labs(title = "Forestry Bog, r2 = 0.02", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_manual(values = c("#d8b365", "#5ab4ac")) + theme(legend.position = "none")
#adonis(x ~ layer, as(sample_data(FB), "data.frame")) #0.01582, 0.101

WS <- prune_samples(sampledata$Bog == "WS" & sampledata$Layer != "U", alldata)
x <- UniFrac(WS, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
scores <- scores(pcoa)
layer <- substr(labels(x), start = 3, stop = 3)
plot.pcoa <- data.frame(scores$sites, layer)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Layer")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
plot3 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Layer)) + geom_point(size=2, alpha = 0.5) + labs(title = "West Sparkling Bog, r2 = 0.12", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_manual(values = c("#d8b365", "#5ab4ac")) + theme(legend.position = "none")
#adonis(x ~ layer, as(sample_data(WS), "data.frame")) #0.10782, 0.001

NS <- prune_samples(sampledata$Bog == "NS" & sampledata$Layer != "U", alldata)
x <- UniFrac(NS, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
scores <- scores(pcoa)
layer <- substr(labels(x), start = 3, stop = 3)
plot.pcoa <- data.frame(scores$sites, layer)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Layer")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
plot4 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Layer)) + geom_point(size=2, alpha = 0.5) + labs(title = "North Sparkling Bog, r2 = 0.19", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_manual(values = c("#d8b365", "#5ab4ac")) + theme(legend.position = "none")
#adonis(x ~ layer, as(sample_data(NS), "data.frame")) #0.18548 p = 0.001

TB <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer != "U", alldata)
x <- UniFrac(TB, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
scores <- scores(pcoa)
layer <- substr(labels(x), start = 3, stop = 3)
plot.pcoa <- data.frame(scores$sites, layer)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Layer")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
plot5 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Layer)) + geom_point(size=2, alpha = 0.5) + labs(title = "Trout Bog, r2 = 0.18", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_manual(values = c("#d8b365", "#5ab4ac")) + theme(legend.position = "none")
#adonis(x ~ layer, as(sample_data(TB), "data.frame")) #0.17675 p = 0.001

SS <- prune_samples(sampledata$Bog == "SS" & sampledata$Layer != "U", alldata)
x <- UniFrac(SS, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
scores <- scores(pcoa)
layer <- substr(labels(x), start = 3, stop = 3)
plot.pcoa <- data.frame(scores$sites, layer)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Layer")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
plot6 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Layer)) + geom_point(size=2, alpha = 0.5) + labs(title = "South Sparkling Bog, r2 = 0.16", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_manual(values = c("#d8b365", "#5ab4ac")) + theme(legend.position = "none")
#adonis(x ~ layer, as(sample_data(SS), "data.frame")) #0.16416 p = 0.001

HK <- prune_samples(sampledata$Bog == "HK" & sampledata$Layer != "U", alldata)
x <- UniFrac(HK, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
scores <- scores(pcoa)
layer <- substr(labels(x), start = 3, stop = 3)
plot.pcoa <- data.frame(scores$sites, layer)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Layer")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
plot7 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Layer)) + geom_point(size=2, alpha = 0.5) + labs(title = "Hell's Kitchen, r2 = 0.40", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_manual(values = c("#d8b365", "#5ab4ac")) + theme(legend.position = "none")
#adonis(x ~ layer, as(sample_data(HK), "data.frame")) #0.40011 p = 0.001

MA <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer != "U", alldata)
x <- UniFrac(MA, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
scores <- scores(pcoa)
layer <- substr(labels(x), start = 3, stop = 3)
plot.pcoa <- data.frame(scores$sites, layer)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Layer")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
plot8 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Layer)) + geom_point(size=2, alpha = 0.5) + labs(title = "Mary Lake, r2 = 0.45", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_manual(values = c("#d8b365", "#5ab4ac")) + theme(legend.position = "none")
#adonis(x ~ layer, as(sample_data(MA), "data.frame")) #0.4505 p = 0.001

layer_splits <- plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, align = "h", nrow = 3)
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots_2017/FigS4.pdf", layer_splits, base_aspect_ratio = 1.5, base_height = 8)

##########
#Figure S6 - alternative colorations of Fig 2

# Separate epilimnion and hypolimnion samples
epi <- prune_samples(sampledata$Layer == "E", alldata)
hypo <- prune_samples(sampledata$Layer == "H", alldata)

# Analyze and plot epilimnion points
#epi <- prune_taxa(taxa_sums(alldata) > 1000, alldata)
x.epi <- UniFrac(epi, weighted = T, normalize = T)
pcoa <- betadisper(x.epi, substr(labels(x.epi), start = 1, stop = 3))
scores <- scores(pcoa)
lakes <- factor(substr(rownames(scores$sites), start = 1, stop = 2), levels = c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA"))
regime <- c()
regime[which(lakes == "CB" | lakes == "FB" | lakes == "WS")] <- "polymictic"
regime[which(lakes == "TB" | lakes == "NS" | lakes == "SS")] <- "dimictic"
regime[which(lakes == "HK" | lakes == "MA")] <- "meromictic"

plot.pcoa <- data.frame(scores$sites, lakes)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Lake")
plot.pcoa$Date <- extract_date(rownames(plot.pcoa))
plot.pcoa$Regime <- regime

# Add in temperature where available
meta2 <- metadata[,c(1,2,3,4)]
meta2$Layer <- substr(meta2$Sample_Name, start = 3, stop = 3)
meta2 <- meta2[which(meta2$Layer == "E"), ]
meta2$Site <- substr(meta2$Sample_Name, start = 1, stop = 2)
meta2$Date <- extract_date(meta2$Sample_Name)

layer <- c()
for(i in 1:dim(meta2)[1]){
  sample <- meta2[i, ]
  if(sample$Site == "CB" | sample$Site == "FB"){
    if(sample$Depth <= 1){
      layer[i] <- "Epi"
    }else if(sample$Depth > 1){
      layer[i] <- "Hypo"
    }
  }else if (sample$Site == "WS" | sample$Site == "NS" | sample$Site == "SS" | sample$Site == "TB"| sample$Site == "HK" | sample$Site == "MA"){
    if(sample$Depth <= 2){
      layer[i] <- "Epi"
    }else if(sample$Depth > 2){
      layer[i] <- "Hypo"
    }
  }
}
meta2$Layer <- layer
meta2$Depth <- NULL
mean_meta2 <- aggregate(x = meta2, by = list(meta2$Layer, meta2$Site, meta2$Date), FUN = "mean")
mean_meta2 <- mean_meta2[which(is.na(mean_meta2$Temperature) == F), ]

tempdata <- c()
for(i in 1:dim(plot.pcoa)[1]){
  sample <- plot.pcoa[i, ]
  hit <- which(mean_meta2$Group.2 == sample$Lake & mean_meta2$Date == sample$Date & mean_meta2$Group.1 == "Epi")
  if(length(hit) > 0){
    tempdata[i] <- mean_meta2$Temperature[hit]
  }else{
    tempdata[i] <- NA
  }
}

plot.pcoa$Temp <- tempdata

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

figS6a <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = as.integer(Date))) + geom_point(size = 1, alpha = 1/2) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_gradientn(colors = c("#ffffcc", "#a1dab4","#41b6c4", "#225ea8"))

figS6b <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Regime)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_blank())

figS6c <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Temp)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_gradientn(colors = c("#2c7bb6", "#abd9e9","#fdae61", "#d7191c")) + theme(legend.title = element_blank())


# Analyze and plot hypolimnion points
x.hypo <- UniFrac(hypo, weighted = T, normalize = T)
pcoa <- betadisper(x.hypo, substr(labels(x.hypo), start = 1, stop = 3))
scores <- scores(pcoa)
years <- factor(substr(rownames(scores$sites), start = 9, stop = 10), levels = c("05", "07", "08", "09"))
lakes <- factor(substr(rownames(scores$sites), start = 1, stop = 2), levels = c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA"))
regime <- c()
regime[which(lakes == "CB" | lakes == "FB" | lakes == "WS")] <- "polymictic"
regime[which(lakes == "TB" | lakes == "NS" | lakes == "SS")] <- "dimictic"
regime[which(lakes == "HK" | lakes == "MA")] <- "meromictic"
plot.pcoa <- data.frame(scores$sites, lakes)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2",  "Lake")
plot.pcoa$Date <- extract_date(rownames(plot.pcoa))
plot.pcoa$Regime <- regime

# Calculate percent variation explained of each axis
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

tempdata <- c()
for(i in 1:dim(plot.pcoa)[1]){
  sample <- plot.pcoa[i, ]
  hit <- which(mean_meta2$Group.2 == sample$Lake & mean_meta2$Date == sample$Date & mean_meta2$Group.1 == "Hypo")
  if(length(hit) > 0){
    tempdata[i] <- mean_meta2$Temperature[hit]
  }else{
    tempdata[i] <- NA
  }
}

plot.pcoa$Temp <- tempdata

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

figS6d <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = as.integer(Date))) + geom_point(size = 1, alpha = 1/2) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_gradientn(colors = c("#ffffcc", "#a1dab4","#41b6c4", "#225ea8"))

figS6e <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Regime)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_blank())

figS6f <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Temp)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_gradientn(colors = c("#2c7bb6", "#abd9e9","#fdae61", "#d7191c")) + theme(legend.title = element_blank())

figS6 <- plot_grid(figS6a,figS6d, figS6b, figS6e, figS6c, figS6f, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"), align = "w")

save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS6.pdf", figS7, base_aspect_ratio = 1.5, base_height = 6)

##########
#Figure S7 - alternative metric PCoA in Fig 2
# Use Bray-Curtis to measure beta diversity between sites
x <- vegdist(t(otu_table), method = "bray")
beta_diversity <- as.matrix(x)

# Compare sites in epilimnion
epis <- c("CBE", "FBE", "WSE", "NSE", "TBE", "SSE", "HKE", "MAE")
epi_comparison <- expand.grid(epis, epis)
epi_results <- c()
for(i in 1:dim(epi_comparison)[1]){
  site1 <- grep(epi_comparison$Var1[i], rownames(beta_diversity))
  site2 <- grep(epi_comparison$Var2[i], colnames(beta_diversity))
  compare <- beta_diversity[site1, site2]
  epi_results[i] <- mean(as.matrix(compare))
}
epi_comparison$Beta_Diversity <- epi_results

p1 <- ggplot(data = epi_comparison, aes(x = Var1, y = Var2, fill = Beta_Diversity)) + geom_tile() + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred", midpoint = 0.7, name = "Bray-Curtis") + labs(x = NULL, y = NULL, title ="Epilimnia") 


hypos <- c("CBH", "FBH", "WSH", "NSH", "TBH", "SSH", "HKH", "MAH")
hypo_comparison <- expand.grid(hypos, hypos)
hypo_results <- c()
for(i in 1:dim(hypo_comparison)[1]){
  site1 <- grep(hypo_comparison$Var1[i], rownames(beta_diversity))
  site2 <- grep(hypo_comparison$Var2[i], colnames(beta_diversity))
  compare <- beta_diversity[site1, site2]
  hypo_results[i] <- mean(as.matrix(compare))
}
hypo_comparison$Beta_Diversity <- hypo_results

p2 <- ggplot(data = hypo_comparison, aes(x = Var1, y = Var2, fill = Beta_Diversity)) + geom_tile() + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred", midpoint = 0.7, name = "Bray-Curtis") + labs(x = NULL, y = NULL, title ="Hypolimnia")
figS7 <- plot_grid(p1, p2, nrow = 1)
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS7.pdf", figS7, base_aspect_ratio = 3, base_height = 4)



##########
#Figure S8 - PCoAs of layers by years not shown in Fig 3
colors <- c("#a6cee3", "#1f78b4", "#33a02c", "#b2df8a")

NSE <- prune_samples(sampledata$Bog == "NS" & sampledata$Layer == "E", alldata)
NSE_year <- factor(substr(sample_names(NSE), start = 9, stop = 10), levels = years)
x <- UniFrac(NSE, weighted = T, normalize = T)
pcoa <- betadisper(x, NSE_year)
scores <- scores(pcoa)
# Locate centroids
NSEcentroids <- scores$centroids
NSEcentroids <- as.data.frame(NSEcentroids)
NSEcentroids$Year <- factor(years[2:4], level = years[2:4])
plot.pcoa <- data.frame(scores$sites, NSE_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
figS8a <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=NSEcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "North Sparkling Epilimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.15, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
adonis(x ~ Year, as(sample_data(NSE), "data.frame"))

TBE <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "E", alldata)
TBE_year <- factor(substr(sample_names(TBE), start = 9, stop = 10), levels = years)
x <- UniFrac(TBE, weighted = T, normalize = T)
pcoa <- betadisper(x, TBE_year)
scores <- scores(pcoa)
# Locate centroids
TBEcentroids <- scores$centroids
TBEcentroids <- as.data.frame(TBEcentroids)
TBEcentroids$Year <- factor(years, level = years)
plot.pcoa <- data.frame(scores$sites, TBE_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
figS8b <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=TBEcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Trout Bog Epilimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.11, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
adonis(x ~ Year, as(sample_data(TBE), "data.frame"))

SSE <- prune_samples(sampledata$Bog == "SS" & sampledata$Layer == "E", alldata)
SSE_year <- factor(substr(sample_names(SSE), start = 9, stop = 10), levels = years)
x <- UniFrac(SSE, weighted = T, normalize = T)
pcoa <- betadisper(x, SSE_year)
scores <- scores(pcoa)
# Locate centroids
SSEcentroids <- scores$centroids
SSEcentroids <- as.data.frame(SSEcentroids)
SSEcentroids$Year <- factor(years[2:4], level = years[2:4])
plot.pcoa <- data.frame(scores$sites, SSE_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
figS8c <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=SSEcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "South Sparkling Epilimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.10, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
adonis(x ~ Year, as(sample_data(SSE), "data.frame"))

MAE <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer == "E", alldata)
MAE_year <- factor(substr(sample_names(MAE), start = 9, stop = 10), levels = years)
x <- UniFrac(MAE, weighted = T, normalize = T)
pcoa <- betadisper(x, MAE_year)
scores <- scores(pcoa)
# Locate centroids
MAEcentroids <- scores$centroids
MAEcentroids <- as.data.frame(MAEcentroids)
MAEcentroids$Year <- factor(years, level = years)
plot.pcoa <- data.frame(scores$sites, MAE_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
figS8d <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=MAEcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Mary Lake Epilimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.12, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
adonis(x ~ Year, as(sample_data(MAE), "data.frame"))

NSH <- prune_samples(sampledata$Bog == "NS" & sampledata$Layer == "H", alldata)
NSH_year <- factor(substr(sample_names(NSH), start = 9, stop = 10), levels = years)
NSH_year[17] = "08"
x <- UniFrac(NSH, weighted = T, normalize = T)
pcoa <- betadisper(x, NSH_year)
scores <- scores(pcoa)
# Locate centroids
NSHcentroids <- scores$centroids
NSHcentroids <- as.data.frame(NSHcentroids)
NSHcentroids$Year <- factor(years[2:4], level = years[2:4])
plot.pcoa <- data.frame(scores$sites, NSH_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")
axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)
figS8e <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=NSHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "North Sparkling Hypolimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.11, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
adonis(x ~ Year, as(sample_data(NSH), "data.frame"))

figS8 <- plot_grid(figS8a, figS8b, figS8c, figS8d, figS8e, align = "h", nrow = 3, labels = c("A", "B", "C", "D", "E"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS8.pdf", figS8, base_aspect_ratio = 1, base_height = 7)

##########
#Figure 9
# Time decay plots

x <- UniFrac(alldata, weighted = T, normalize = T)
beta_diversity <- as.matrix(x)

site <- c("CBE", "FBE", "WSE", "NSE", "TBE", "SSE", "HKE", "MAE", "CBH", "FBH", "WSH", "NSH", "TBH", "SSH", "HKH", "MAH")

# For each site, make plots of beta diversity of tn:t0 and tn:tn-1

for(i in 1:length(site)){
  rows <- grep(site[i], rownames(beta_diversity))
  cols <- grep(site[i], colnames(beta_diversity))
  lake <- beta_diversity[rows,cols]
  datekey <- extract_date(colnames(lake))
  
  #Compare to time 0
  start <- lake[which(datekey == min(datekey)), ]
  if(length(which(datekey == min(datekey))) > 1){
    start <- start[1, ]
  }
  datediff <- datekey - min(datekey)
  plot_decay <- data.frame(as.numeric(start), as.numeric(datediff))
  colnames(plot_decay) <- c("UniFrac", "Days")
  plot_decay <- plot_decay[order(plot_decay$Days), ]
  object1 <- ggplot(data = plot_decay, aes(x = Days, y = UniFrac)) + geom_point() + geom_path() + labs(x = "Days since 1st Timepoint", y = "Weighted UniFrac", title = paste(site[i], "Time Decay"))
  print(object1)
  
  #Compare to previous timepoint
  lake <- lake[order(datekey), order(datekey)]
  lake <- lake[grep(".R2", rownames(lake), invert = T), grep(".R2", colnames(lake), invert = T)]
  t1 <- colnames(lake)[1:length(colnames(lake))-1]
  t2 <- colnames(lake)[2:length(colnames(lake))]
  set <- c()
  for(k in 1:length(t1)){
    set[k] <- lake[grep(t1[k], rownames(lake)), grep(t2[k], colnames(lake))]
  }
  plot_turnover <- data.frame(t1, t2, set)
  colnames(plot_turnover) <- c("Sample1", "Sample2", "Diversity")
  plot_turnover$Date1 <- extract_date(plot_turnover$Sample1)
  plot_turnover$Date2 <- extract_date(plot_turnover$Sample2)
  plot_turnover$Days <- as.numeric(plot_turnover$Date2 - plot_turnover$Date1)
  plot_turnover$NormDiversity <- plot_turnover$Diversity/plot_turnover$Days
  
  object2 <- ggplot(data = plot_turnover, aes(x = Date2, y = NormDiversity)) + geom_point() + geom_path() + labs(x = c("Date"), y = c("Turnover"), title = paste(site[i], "Turnover"))
  print(object2)
}

##########
#FigS10 - zscore normalized OTUs over multiple years

annual_trends <- function(lake, otu){
  bog <- bog_subset(lake, otu_table)
  year1 <- year_subset("05", bog)
  year2 <- year_subset("07", bog)
  year3 <- year_subset("08", bog)
  year4 <- year_subset("09", bog)
  
  if(dim(year1)[2] > 0){
    year1 <- zscore(year1)
    year2 <- zscore(year2)
    year3 <- zscore(year3)
    year4 <- zscore(year4)
    
    ztable <- cbind(year1, year2, year3, year4)
  }else if(dim(year1)[2] == 0 & dim(year3)[2] > 0){
    year2 <- zscore(year2)
    year3 <- zscore(year3)
    year4 <- zscore(year4)
    
    ztable <- cbind(year2, year3, year4)
  }else if(dim(year1)[2] == 0 & dim(year3)[2] == 0 & dim(year4)[2] > 0){
    year2 <- zscore(year2)
    year4 <- zscore(year4)
    
    ztable <- cbind(year2, year4)
  }else{
    ztable <- zscore(year2)
  }
  ztable <- melt(ztable)
  ztable$Year <- substr(ztable$Var2, start = 9, stop = 10)
  ztable$Day <- as.numeric(format(extract_date(ztable$Var2), format = "%j"))
  
  plot <- ggplot(data = ztable[which(ztable$Var1 == otu), ], aes(x = Day, y = value, group = Year, color = Year)) + geom_point() + geom_line() + theme_bw() + labs(title = paste(lake, otu, sep = ": "))
  return(plot)
}

plot1 <- annual_trends("TBE", "Otu0002")
plot2 <- annual_trends("TBH", "Otu0002")
plot3 <- annual_trends("TBH", "Otu0050")
plot4 <- annual_trends("MAH", "Otu0050")

figS10 <- plot_grid(plot1, plot2, plot3, plot4, align = "h", nrow = 2)
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS10.pdf", figS10, base_aspect_ratio = 2.5, base_height = 4)

##########
#Figure S11 - rarefaction curves

rarefaction <- function(table){
  groups <- c()
  totals <- c()
  table <- table[, sample(ncol(table))]
  for(i in 1:dim(table)[2]){
    groups <- append(groups, rownames(table)[which(table[,i] > 0)], length(groups))
    totals[i] <- length(unique(groups))
  }
  return(totals)
}

all.rf <- rarefaction(otu_table)
all.rf <- data.frame(seq(1:length(all.rf)), all.rf)
colnames(all.rf) <- c("Index", "OTUs")
figS11a <- ggplot(all.rf, aes(x = Index, y = OTUs)) + geom_line(size = 1) + labs(title = "All Data", x = "Number of Samples", y = "Number of Taxa")

epi.rf <- rarefaction(bog_subset("..E", otu_table))
epi.rf <- data.frame(seq(1:length(epi.rf)), epi.rf)
colnames(epi.rf) <- c("Index", "OTUs")
figS11b <- ggplot(epi.rf, aes(x = Index, y = OTUs)) + geom_line(size = 1) + labs(title = "Epilimnion Data", x = "Number of Samples", y = "Number of Taxa")

hypo.rf <- rarefaction(bog_subset("..H", otu_table))
hypo.rf <- data.frame(seq(1:length(hypo.rf)), hypo.rf)
colnames(hypo.rf) <- c("Index", "OTUs")
figS11c <- ggplot(hypo.rf, aes(x = Index, y = OTUs)) + geom_line(size = 1) + labs(title = "Hypo Data", x = "Number of Samples", y = "Number of Taxa")

figS11 <- plot_grid(figS11a, figS11b, figS11c, align = "h", nrow = 1, labels = c("A", "B", "C"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS11.pdf", figS11, base_aspect_ratio = 3, base_height = 3)


### Figure S8 - Consistent lineage traits by year
lineage_table <- combine_otus("Lineage", otu_table, taxonomy)

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
  
  z <- ggplot(data = to.plot, aes(x = variance, y = abundance, fill = persistance)) + geom_point() + theme_bw() + labs(title = paste(year[i]), x = "Variability (CV)", y = "Mean Abundance when Present") + geom_label_repel(aes(label = lineage)) + scale_fill_gradient(low = "white", high = "lightgreen")  + theme(legend.title = element_blank())
  
  assign(paste("p", i, sep = ""), z)
}

figS12_legend <- get_legend(p2)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

figS12 <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = c("A", "B", "C", "D"))
figS12 <- plot_grid(figS12, figS12_legend, nrow = 1, rel_widths = c(1, 0.1))

save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/FigS12.pdf", figS12, base_aspect_ratio = 1.5, base_height = 6)

# ##########
# #Indicator analysis
# 
# # Make tables at each phylogenetic level
# tribe_table <- combine_otus("Tribe", otu_table, taxonomy)
# clade_table <- combine_otus("Clade", otu_table, taxonomy)
# lineage_table <- combine_otus("Lineage", otu_table, taxonomy)
# order_table <- combine_otus("Order", otu_table, taxonomy)
# class_table <- combine_otus("Class", otu_table, taxonomy)
# phylum_table <- combine_otus("Phylum", otu_table, taxonomy)
# 
# # Change OTU number designations to full taxonomic assignment
# named_otu_table <- otu_table
# fullnames <- c()
# for(i in 1:dim(taxonomy)[1]){
#   fullnames[i] <- paste(taxonomy[i, ], collapse = ";")
# }
# fullnames <- make.unique(fullnames)
# rownames(named_otu_table) <- fullnames
# 
# # Combine tables at each level into one giant table
# full_table <- rbind(named_otu_table, tribe_table, clade_table, lineage_table, order_table, class_table, phylum_table)
# 
# # Remove groups that are unclassified at any level
# classified <- grep("unclassified", rownames(full_table))
# classified1 <- grep("__$", rownames(full_table))
# input_table <- full_table[-c(classified, classified1),]
# 
# # Keep only the top quantile in abundance
# threshold <- 500
# input_table <- input_table[which(rowSums(input_table) >= threshold),]
# 
# # Format for input into mulitplatt() function
# input_table <- t(input_table)
# input_table <- as.data.frame(input_table)
# 
# # Group by mixing regime
# lakeid <- c("CB", "FB", "WS", "NS", "TB", "SS", "HK","MA")
# lakes <- substr(rownames(input_table), start = 1, stop = 2)
# 
# mixing_groups <- c()
# mixing_groups[which(lakes == "CB" | lakes == "FB" | lakes == "WS")] <- 1
# mixing_groups[which(lakes == "TB"| lakes == "NS"| lakes == "SS")] <- 2
# mixing_groups[which(lakes == "MA" | lakes == "HK")] <- 3
# 
# # Run indicator taxa analysis
# clades_by_mixing <- multipatt(x = input_table, cluster = mixing_groups, func = "r.g", control = how(nperm = 9999))
# 
# # Save data
# write.csv(clades_by_mixing$sign, file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/indicators_by_mixing_2017-01-17.csv")




