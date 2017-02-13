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
#library(indicspecies)

data(metadata)
#data(otu_table)
otu_table <- read.csv("C:/Users/Alex/Documents/deblurred_nonrarified_otutable.csv", header = T, row.names = 1)
#data(taxonomy)
taxonomy <- read.csv("C:/Users/Alex/Documents/nonrarified_clean_taxonomy.csv", header = T, row.names = 1, colClasses = c("character"))
seqs2 <- read.dna("C:/Users/Alex/Documents/bog_repseqs_nonrarified_aligned.good.filter.fasta", format = "fasta", as.matrix = T)
otu_table <- otu_table[match(dimnames(seqs2)[[1]], rownames(otu_table)), ]
taxonomy <- taxonomy[match(rownames(otu_table), rownames(taxonomy)), ]
d <- dist.dna(seqs2, model = "raw")
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
fig1a <- ggplot(data = epi.data, aes(y = epi.obs, x = epi.lakes, fill = epi.lakes)) + geom_boxplot() + labs(y = "Observed Richness", x = NULL) + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_blank()) + labs(title = "Epilimnion") + theme(legend.position = "none")
fig1a <- fig1a + annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 550, ymax = 537.5, fill = "#a6cee3", alpha = 0.5) #CB1
fig1a <- fig1a + annotate("rect", xmin = 4.5, xmax = 8.5, ymin = 550, ymax = 537.5, fill = "#a6cee3", alpha = 0.5) #CB2
fig1a <- fig1a + annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 525, ymax = 537.5, fill = "#1f78b4", alpha = 0.5) #FB1
fig1a <- fig1a + annotate("rect", xmin = 2.5, xmax = 8.5, ymin = 525, ymax = 537.5, fill = "#1f78b4", alpha = 0.5) #FB2
fig1a <- fig1a + annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 512.5, ymax = 525, fill = "#b2df8a", alpha = 0.5) #WS1
fig1a <- fig1a + annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 512.5, ymax = 525, fill = "#b2df8a", alpha = 0.5) #WS2
fig1a <- fig1a + annotate("rect", xmin = 0.5, xmax = 3.5, ymin = 500, ymax = 512.5, fill = "#33a02c", alpha = 0.5) #NS1
fig1a <- fig1a + annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 500, ymax = 512.5, fill = "#33a02c", alpha = 0.5) #NS2
fig1a <- fig1a + annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 487.5, ymax = 500, fill = "#fb9a99", alpha = 0.5) #TB1
fig1a <- fig1a + annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 487.5, ymax = 500, fill = "#fb9a99", alpha = 0.5) #TB2
fig1a <- fig1a + annotate("rect", xmin = 0.5, xmax = 5.5, ymin = 475, ymax = 487.5, fill = "#e31a1c", alpha = 0.5) #SS1
fig1a <- fig1a + annotate("rect", xmin = 7.5, xmax = 8.5, ymin = 475, ymax = 487.5, fill = "#e31a1c", alpha = 0.5) #SS2
fig1a <- fig1a + annotate("rect", xmin = 0.5, xmax = 5.5, ymin = 462.5, ymax = 475, fill = "#fdbf6f", alpha = 0.5) #HK1
fig1a <- fig1a + annotate("rect", xmin = 7.5, xmax = 8.5, ymin = 462.5, ymax = 475, fill = "#fdbf6f", alpha = 0.5) #HK2
fig1a <- fig1a + annotate("rect", xmin = 0.5, xmax = 7.5, ymin = 450, ymax = 462.5, fill = "#ff7f00", alpha = 0.5) #MA


# plot 1B
fig1b <- ggplot(data = hypo.data, aes(y = hypo.obs, x = hypo.lakes, fill = hypo.lakes)) + geom_boxplot() + labs(y="Observed Richness", x = NULL) + scale_fill_brewer(palette = "Paired") + theme(legend.position = "none") + labs(title = "Hypolimnion")
fig1b <- fig1b + annotate("rect", xmin = 1.5, xmax = 8.5, ymin = 650, ymax = 637.5, fill = "#a6cee3", alpha = 0.5) #CB1
fig1b <- fig1b + annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 625, ymax = 637.5, fill = "#1f78b4", alpha = 0.5) #FB1
fig1b <- fig1b + annotate("rect", xmin = 2.5, xmax = 8.5, ymin = 625, ymax = 637.5, fill = "#1f78b4", alpha = 0.5) #FB2
fig1b <- fig1b + annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 612.5, ymax = 625, fill = "#b2df8a", alpha = 0.5) #WS1
fig1b <- fig1b + annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 612.5, ymax = 625, fill = "#b2df8a", alpha = 0.5) #WS2
fig1b <- fig1b + annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 600, ymax = 612.5, fill = "#33a02c", alpha = 0.5) #NS1
fig1b <- fig1b + annotate("rect", xmin = 6.5, xmax = 8.5, ymin = 600, ymax = 612.5, fill = "#33a02c", alpha = 0.5) #NS2
fig1b <- fig1b + annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 587.5, ymax = 600, fill = "#fb9a99", alpha = 0.5) #TB1
fig1b <- fig1b + annotate("rect", xmin = 6.5, xmax = 8.5, ymin = 587.5, ymax = 600, fill = "#fb9a99", alpha = 0.5) #TB2
fig1b <- fig1b + annotate("rect", xmin = 0.5, xmax = 5.5, ymin = 575, ymax = 587.5, fill = "#e31a1c", alpha = 0.5) #SS1
fig1b <- fig1b + annotate("rect", xmin = 6.5, xmax = 8.5, ymin = 575, ymax = 587.5, fill = "#e31a1c", alpha = 0.5) #SS2
fig1b <- fig1b + annotate("rect", xmin = 0.5, xmax = 6.5, ymin = 562.5, ymax =575, fill = "#fdbf6f", alpha = 0.5) #HK1
fig1b <- fig1b + annotate("rect", xmin = 7.5, xmax = 8.5, ymin = 562.5, ymax = 575, fill = "#fdbf6f", alpha = 0.5) #HK2
fig1b <- fig1b + annotate("rect", xmin = 0.5, xmax = 7.5, ymin = 550, ymax = 562.5, fill = "#ff7f00", alpha = 0.5) #MA

fig1 <- plot_grid(fig1a, fig1b, nrow = 2, labels = c("A", "B"), align = "w")
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/Fig1nonrarified.pdf", fig1, base_aspect_ratio = 1.3, base_height = 4.75)


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

fig2a <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + geom_path(data = df_ell, aes(x = PCoA1, y = PCoA2, colour = Lake), size = 1, linetype = 1) + theme(legend.title = element_blank())

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

fig2b <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + geom_path(data = df_ell, aes(x = PCoA1, y = PCoA2, colour = Lake), size = 1, linetype = 1) + theme(legend.position = "none")


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
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/Fig2nonrarified.pdf", fig2, base_aspect_ratio = 1.2, base_height = 6)

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

TBH_year <- TBH_year[which(is.na(TBH_year) == F)]
plot.pcoa <- data.frame(scores$sites, TBH_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

fig3a <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + geom_point(data=TBHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Trout Bog", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.36, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors) + theme(plot.subtitle = element_text(hjust = 0.5))
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

fig3b <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + geom_point(data=SSHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "South Sparkling Bog", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.20, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors[2:4])  + theme(legend.position = "none", plot.subtitle = element_text(hjust = 0.5))

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

fig3c <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=MAHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Mary Lake", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.09, p = 0.002") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)

# Calculate PERMANOVA - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
adonis(x ~ Year, as(sample_data(MAH), "data.frame"))
# r2 0.09169
# Pr 0.002

#Panel 2

# calculate unifrac, then run the betadisper() function on different groupings of the data and boxplot their dispersion
x <- UniFrac(alldata, weighted = T, normalize = T)
#Rerun with this x to remove outliers
#x <- UniFrac(prune_samples(rownames(sampledata) != rownames(plot.boxes)[which(plot.boxes$Distance > 0.40)], alldata), weighted = T, normalize = T)
#layer_group <- betadisper(x, substr(labels(x), start = 3, stop = 3))
lakelayer_group <- betadisper(x, substr(labels(x), start = 1, stop = 3))
plot.boxes <- data.frame(lakelayer_group$group, lakelayer_group$distances)
colnames(plot.boxes) <- c("Group", "Distance")
plot.boxes <- plot.boxes[which(plot.boxes$Group != "CBU" & plot.boxes$Group != "NSU"), ]
plot.boxes$Group <- factor(plot.boxes$Group, levels = c("MAE", "MAH",  "HKE", "HKH", "SSE", "SSH", "TBE", "TBH",  "NSE", "NSH",  "WSE", "WSH", "FBE", "FBH", "CBE", "CBH"))

lakelayerkey <- c("MAE", "MAH", "HKE", "HKH", "SSE", "SSH", "TBE", "TBH",  "NSE", "NSH",  "WSE", "WSH", "FBE", "FBH", "CBE", "CBH")
x1 <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5)
pal <- rev(c("#a6cee3", "#a6cee3", "#1f78b4", "#1f78b4", "#b2df8a", "#b2df8a", "#33a02c", "#33a02c", "#fb9a99", "#fb9a99", "#e31a1c", "#e31a1c", "#fdbf6f", "#fdbf6f", "#ff7f00", "#ff7f00"))
background <- data.frame(lakelayerkey, x1, pal)

fig3d <- ggplot() + geom_boxplot(data = plot.boxes, aes(x = Group, y = Distance)) + coord_flip() + annotate("rect", xmin = x1, xmax = x1 + 1, ymin = -Inf, ymax = Inf, fill = pal) + geom_boxplot(data = plot.boxes, aes(x = Group, y = Distance)) + labs(y = "Dispersion", x = NULL)

fig3_panel1 <- plot_grid(fig3a, fig3b, fig3c, nrow = 3, labels = c("A", "B", "C"))
fig3_panel1 <- plot_grid(fig3_panel1, fig3_legend, nrow = 1, rel_widths = c(1, 0.1))
fig3 <- plot_grid(fig3_panel1, fig3d, nrow = 1, rel_widths = c(1, 1), labels = c("", "D"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/Fig3nonrarified.pdf", fig3, base_aspect_ratio = 1.5, base_height = 6)


# Do layers in the same lake have significantly different dispersion?
permutest(lakelayer_group, pairwise = T)
# Mulitple p-values by 8 (my number of comparisons) for a Bonferroni correction
# Only TB and SS significant


###############
#Figure 4
#Venn Diagrams

#Part 1 - core analysis, in text

# Is there a natural breaking point?
# FWcore <- rownames(otu_table)[which(rowSums(otu_table > 0) == 1387)]
# FWcore <- taxonomy[match(FWcore, rownames(taxonomy)), ]
# 
# persistant_levels <- c()
# for(i in 1:100){
#   persistant_levels[i] <- length(which(rowSums(otu_table > 0) >= 1387*(i/100)))
# }
# barplot(persistant_levels)
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

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/Fig4anonrarified.pdf", width = 2, height = 2)
grid.draw(v1)
dev.off()

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/Fig4bnonrarified.pdf", width = 2, height = 2)
grid.draw(v2)
dev.off()

##########
#Figure 5 - the richness over time trace

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

fig5a <- ggplot() + geom_line(data = TB_richness, aes(x = date, y = richness), size = 1) + labs(title = "Trout Bog", x = NULL, y = "Observed Richness") + geom_point(data = TB_richness[match(TBHmixes, TB_richness$date), ], aes(x = date, y = richness), size = 2, colour = "red")


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

fig5b <- ggplot() + geom_line(data = NS_richness, aes(x = date, y = richness), size = 1.2) + labs(title = "North Sparkling Bog", x = NULL, y = "Observed Richness")  + geom_point(data = NS_richness[match(NSHmixes, NS_richness$date), ], aes(x = date, y = richness), size = 2, colour = "red")
fig5 <- plot_grid(fig5a, fig5b, align = "h", nrow = 2, labels = c("A", "B"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/Fig5nonrarified.pdf", fig5, base_aspect_ratio = 1.5, base_height = 5)

# What were the most abundant taxa during both mixes?

TBH <- bog_subset("TBH", otu_table)
TBHdates <- extract_date(colnames(TBH))
TBH_mixed_samples <- TBH[, match(TBHmixes, TBHdates)]
topTBH <- sort(rowSums(TBH_mixed_samples), decreasing = T)[1:10]
topTBH_tax <- taxonomy[match(names(topTBH), rownames(taxonomy)), ]

NSH <- bog_subset("NSH", otu_table)
NSHdates <- extract_date(colnames(NSH))
NSH_mixed_samples <- NSH[, match(NSHmixes, NSHdates)[c(1:4, 6)]]
topNSH <- sort(rowSums(NSH_mixed_samples), decreasing = T)[1:10]
topNSH_tax <- taxonomy[match(names(topNSH), rownames(taxonomy)), ]


##########
#Supplemental figures
##########
#Figure S1
#Phylum rank abundance
phylum_table <- combine_otus("Phylum", otu_table, taxonomy)
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

figs1 <- ggplot(data=phylum_totals, aes(x = phyla.names, y = totals)) + geom_bar(stat = "identity", fill = "grey", colour="black") + labs(x = NULL, y = "Log of Observed Reads") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, colour = "black"), axis.title = element_text(size = 15, vjust=1.2), plot.title = element_text(size = 20), axis.text.y = element_text(colour="black")) + scale_y_log10(expand = c(0,0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/FigS1nonrarified.pdf", figs1, base_aspect_ratio = 2, base_height = 5)

###########
#Fig S2 - Run one plot with all the data
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

figS2 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake, shape = Layer)) + geom_point(size = 3, alpha = 1/2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), panel.border = element_rect(colour = "black", fill=NA, size=1)) + labs(title = "All Data", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_shape_manual(values=c(21, 22, 23, 24))  + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired")
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/FigS2nonrarified.pdf", figS2, base_aspect_ratio = 2.5, base_height = 4)

##########
#Figure S3 - Run the UniFrac by lake by layer
CB <- prune_samples(sampledata$Bog == "CB" & sampledata$Layer != "U", alldata)
x <- UniFrac(CB, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
scores <- scores(pcoa)
layer <- substr(labels(x), start = 3, stop = 3)
plot.pcoa <- data.frame(scores$sites, layer)
plot.pcoa <- plot.pcoa[2:dim(plot.pcoa)[1], ]
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
plot.pcoa$Layer[186] <- "H"
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
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/FigS3nonrarified.pdf", layer_splits, base_aspect_ratio = 1.5, base_height = 8)

##########
#Figure S4 - PCoAs of layers by years not shown in Fig 3

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
figS4a <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=NSEcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "North Sparkling Epilimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.15, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
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
figS4b <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=TBEcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Trout Bog Epilimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.11, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
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
figS4c <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=SSEcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "South Sparkling Epilimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.10, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
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
figS4d <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=MAEcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Mary Lake Epilimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.12, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
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
figS4e <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none", plot.subtitle = element_text(hjust = 0.5)) + geom_point(data=NSHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "North Sparkling Hypolimnion", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"), subtitle = "r2 = 0.11, p = 0.001") + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
adonis(x ~ Year, as(sample_data(NSH), "data.frame"))

figS4 <- plot_grid(figS4a, figS4b, figS4c, figS4d, figS4e, align = "h", nrow = 3, labels = c("A", "B", "C", "D", "E"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/FigS4nonrarified.pdf", figS4, base_aspect_ratio = 1, base_height = 7)

##########
#FigS5 - zscore normalized OTUs over multiple years

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

figS5 <- plot_grid(plot1, plot2, plot3, plot4, align = "h", nrow = 2)
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/FigS5nonrarified.pdf", figS5, base_aspect_ratio = 2.5, base_height = 4)

##########
#Figure S6 - rarefaction curves

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
figS6a <- ggplot(all.rf, aes(x = Index, y = OTUs)) + geom_line(size = 1) + labs(title = "All Data", x = "Number of Samples", y = "Number of Taxa")

epi.rf <- rarefaction(bog_subset("..E", otu_table))
epi.rf <- data.frame(seq(1:length(epi.rf)), epi.rf)
colnames(epi.rf) <- c("Index", "OTUs")
figS6b <- ggplot(epi.rf, aes(x = Index, y = OTUs)) + geom_line(size = 1) + labs(title = "Epilimnion Data", x = "Number of Samples", y = "Number of Taxa")

hypo.rf <- rarefaction(bog_subset("..H", otu_table))
hypo.rf <- data.frame(seq(1:length(hypo.rf)), hypo.rf)
colnames(hypo.rf) <- c("Index", "OTUs")
figS6c <- ggplot(hypo.rf, aes(x = Index, y = OTUs)) + geom_line(size = 1) + labs(title = "Hypo Data", x = "Number of Samples", y = "Number of Taxa")

figS6 <- plot_grid(figS6a, figS6b, figS6c, align = "h", nrow = 1, labels = c("A", "B", "C"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/FigS6nonrarified.pdf", figS6, base_aspect_ratio = 3, base_height = 3)

# # Bonus: nonrandom rarefaction
# rarefaction <- function(table){
#   groups <- c()
#   totals <- c()
#   #table <- table[, sample(ncol(table))]
#   for(i in 1:dim(table)[2]){
#     groups <- append(groups, rownames(table)[which(table[,i] > 0)], length(groups))
#     totals[i] <- length(unique(groups))
#   }
#   return(totals)
# }
# 
# all.rf <- rarefaction(otu_table)
# all.rf <- data.frame(seq(1:length(all.rf)), all.rf)
# colnames(all.rf) <- c("Index", "OTUs")
# bonus <- ggplot(all.rf, aes(x = Index, y = OTUs)) + geom_line(size = 1) + labs(title = "All Data", x = "Number of Samples", y = "Number of Taxa")
# save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/nonrarified_test/bonusnonrarified.pdf", bonus, base_aspect_ratio = 1.2, base_height = 3)


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




