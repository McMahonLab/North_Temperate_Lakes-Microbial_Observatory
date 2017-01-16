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

data(metadata)
data(otu_table)
data(taxonomy)
seqs <- read.dna("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_repseqs_07Jul15.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

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

#pdf(file = paste(path2repo, "epi_boxplot.pdf", sep = ""), width = 3.3125, height = 3)
fig1a <- ggplot(data = epi.data, aes(y = epi.obs, x = epi.lakes, fill = epi.lakes)) + geom_boxplot() + labs(y = "Observed Richness", x = NULL) + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_blank()) + labs(title = "Epilimnion")
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

legend <- get_legend(fig1a)
fig1a <- fig1a + theme(legend.position = "none")

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
fig1 <- plot_grid(fig1, legend, nrow = 1, rel_widths = c(3, 0.3))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots_2017/Fig1.pdf", fig1, base_aspect_ratio = 1.3, base_height = 4.75)


# Check significance using the Wilcoxon Rank Test on medians (not output as pdf, indicated as symbols in Illustrator)
pairwise.wilcox.test(epi.data$epi.obs, epi.data$epi.lakes, p.adjust.method = "bonferroni")

pairwise.wilcox.test(hypo.data$hypo.obs, hypo.data$hypo.lakes, p.adjust.method = "bonferroni")

##########
#Figure 2

# Set up phyloseq object to run UniFrac on
OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(taxonomy))
sampledata <- sample_data(data.frame(Bog = substr(colnames(otu_table), start = 1, stop = 2), Layer = substr(colnames(otu_table), start = 3, stop = 3), Year = substr(colnames(otu_table), start = 9, stop = 10), row.names = colnames(otu_table), stringsAsfactors = F)) 
alldata <- phyloseq(OTU, TAX, sampledata, bogtree)

#Run one plot with all the data
x <- UniFrac(alldata, weighted = T, normalize = T)
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

ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake, shape = Layer)) + geom_point(size = 3, alpha = 1/2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), panel.border = element_rect(colour = "black", fill=NA, size=1)) + labs(title = "All Data", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_shape_manual(values=c(21, 22, 23, 24))  + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired")

#Run the UniFrac by lake by layer, but don't plot - just get the r^2 from ANOSIM
CB <- prune_samples(sampledata$Bog == "CB" & sampledata$Layer != "U", alldata)
x <- UniFrac(CB, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
layer <- substr(labels(x), start = 3, stop = 3)
adonis(x ~ layer, as(sample_data(CB), "data.frame")) #0.03235 p = 0.003

FB <- prune_samples(sampledata$Bog == "FB", alldata)
x <- UniFrac(FB, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
layer <- substr(labels(x), start = 3, stop = 3)
adonis(x ~ layer, as(sample_data(FB), "data.frame")) #0.01592 p = 0.107

WS <- prune_samples(sampledata$Bog == "WS", alldata)
x <- UniFrac(WS, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
layer <- substr(labels(x), start = 3, stop = 3)
adonis(x ~ layer, as(sample_data(WS), "data.frame")) #0.10988 p < 0.001

NS <- prune_samples(sampledata$Bog == "NS" & sampledata$Layer != "U", alldata)
x <- UniFrac(NS, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
layer <- substr(labels(x), start = 3, stop = 3)
adonis(x ~ layer, as(sample_data(NS), "data.frame")) #0.1782, p < 0.001

TB <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer != "U", alldata)
x <- UniFrac(TB, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
layer <- substr(labels(x), start = 3, stop = 3)
adonis(x ~ layer, as(sample_data(TB), "data.frame")) #0.14581 p < 0.001

SS <- prune_samples(sampledata$Bog == "SS" & sampledata$Layer != "U", alldata)
x <- UniFrac(SS, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
layer <- substr(labels(x), start = 3, stop = 3)
adonis(x ~ layer, as(sample_data(SS), "data.frame")) #0.16489 p < 0.001

HK <- prune_samples(sampledata$Bog == "HK" & sampledata$Layer != "U", alldata)
x <- UniFrac(HK, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
layer <- substr(labels(x), start = 3, stop = 3)
adonis(x ~ layer, as(sample_data(HK), "data.frame")) #0.42182 p < 0.001

MA <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer != "U", alldata)
x <- UniFrac(MA, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 3, stop = 3))
layer <- substr(labels(x), start = 3, stop = 3)
adonis(x ~ layer, as(sample_data(MA), "data.frame")) #0.5 p < 0.001

lakekey <- c("CB", "FB", "WS", "NS", "TB", "SS", "HK", "MA")
lakekey <- factor(lakekey, levels = lakekey)
diffs <- c(0.03235, 0.01592, 0.10988, 0.1782, 0.14581, 0.16489, 0.42182, 0.5)
layer_barplot <- data.frame(lakekey, diffs)
ggplot(data = layer_barplot, aes(x = lakekey, y = diffs)) + geom_bar(stat = "identity")

# Separate epilimnion and hypolimnion samples
epi <- prune_samples(sampledata$Layer == "E", alldata)
hypo <- prune_samples(sampledata$Layer == "H", alldata)

# Analyze and plot epilimnion points
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
# Calculate PERMADISP by mixing regime
adonis(x ~ regime, as(cbind(sample_data(hypo), regime), "data.frame"))
# r^2 = 0.21907, p = 0.001

legend <- get_legend(fig2a)
fig2a <- fig2a + theme(legend.position = "none")

fig2 <- plot_grid(fig2a, fig2b, nrow = 2, labels = c("A", "B"), align = "w")
fig2 <- plot_grid(fig2, legend, nrow = 1, rel_widths = c(1, 0.1))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots_2017/Fig2.pdf", fig2, base_aspect_ratio = 1.2, base_height = 6)

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

fig3a <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + geom_point(data=TBHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Trout Bog", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)
fig3_legend <- get_legend(fig3a)
fig3a <- fig3a + theme(legend.position = "none")

# Calculate PERMADISP - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
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

fig3b <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + geom_point(data=SSHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "South Sparkling Bog", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors[2:4])  + theme(legend.position = "none")

# Calculate PERMADISP - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
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

fig3c <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=1) + theme(legend.position="none") + geom_point(data=MAHcentroids, aes(x = PCoA1, y = PCoA2, color = Year), size = 3, shape = 3, color = "black") + labs(title = "Mary Lake", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) + scale_color_manual(values = colors)

# Calculate PERMADISP - I'm interested in Pr (p-value) and R2 (amount of variance explained by year)
adonis(x ~ Year, as(sample_data(MAH), "data.frame"))
# r2 0.09169
# Pr 0.002

#Panel 2

# calculate unifrac, then run the betadisper() function on different groupings of the data and boxplot their dispersion
x <- UniFrac(alldata, weighted = T, normalize = T)
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
save_plot("C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots_2017/Fig3.pdf", fig3, base_aspect_ratio = 1.5, base_height = 6)


# Do layers in the same lake have significantly different dispersion?
pairwise.wilcox.test(lakelayer_group$distances, lakelayer_group$group, p.adjust.method = "bonferroni")
# Not different: CB, FB, HK, NS, WS
# Different: MA, SS, TB
# Conclusion: dimictic/meromictic hypolimnia are less variable than epilimnia

# split by year and repeat (but don't save plots)


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
v1 <- draw.triple.venn(area1 = length(which(poly.epi == T)), area2 = length(which(di.epi == T)), area3 = length(which(mero.epi == T)), n12 = length(which(poly.epi == T & di.epi == T)), n13 = length(which(poly.epi == T & mero.epi == T)), n23 = length(which(mero.epi == T & di.epi == T)), n123 = length(which(poly.epi == T & di.epi == T & mero.epi == T)), category = c("Polymictic", "Dimictic", "Meromictic"), main = c("Epilimnion"), fill = c("#a6cee3", "#33a02c", "#fdbf6f"), alpha = 0.4, cex = 1, cat.cex = 0.45)

grid.newpage()
v2 <- draw.triple.venn(area1 = length(which(poly.hypo == T)), area2 = length(which(di.hypo == T)), area3 = length(which(mero.hypo == T)), n12 = length(which(poly.hypo == T & di.hypo == T)), n13 = length(which(poly.hypo == T & mero.hypo == T)), n23 = length(which(mero.hypo == T & di.hypo == T)), n123 = length(which(poly.hypo == T & di.hypo == T & mero.hypo == T)), category = c("Polymictic", "Dimictic", "Meromictic"), fill = c("#a6cee3", "#33a02c", "#fdbf6f"), alpha = 0.4, cex = 1, cat.cex = 0.45)

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots_2017/Fig4a.pdf", width = 2, height = 2)
grid.draw(v1)
dev.off()

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots_2017/Fig4b.pdf", width = 2, height = 2)
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
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots_2017/Fig5.pdf", fig5, base_aspect_ratio = 1.5, base_height = 5)

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
