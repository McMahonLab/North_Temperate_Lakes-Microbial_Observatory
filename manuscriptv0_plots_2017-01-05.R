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

ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake, fill = Lake, shape = Layer)) + geom_point(size = 3, alpha = 1/2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), panel.border = element_rect(colour = "black", fill=NA, size=1)) + labs(title = "All Data", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_shape_manual(values=c(21, 22, 23, 24))  + scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")

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
