#Analyses in response to reviewers

library(OTUtable)
library(reshape2)
library(ggplot2)
library(ape)
library(phyloseq)
library(vegan)
# library(VennDiagram)
# library(exactRankTests)
library(cowplot)
# library(raster)
# library(ggrepel)
# library(indicspecies)

data(metadata)
data(otu_table)
data(taxonomy)
seqs <- read.dna("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_repseqs_07Jul15.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(taxonomy))
sampledata <- sample_data(data.frame(Bog = substr(colnames(otu_table), start = 1, stop = 2), Layer = substr(colnames(otu_table), start = 3, stop = 3), Year = substr(colnames(otu_table), start = 9, stop = 10), row.names = colnames(otu_table), stringsAsfactors = F)) 
alldata <- phyloseq(OTU, TAX, sampledata, bogtree)


#### Alternative metrics to PCoA

# Use weighted UniFrac to measure beta diversity between sites
x <- UniFrac(alldata, weighted = T, normalize = T)

beta_diversity <- as.matrix(x)
write.csv(beta_diversity, file = "pairwase_beta.csv")

beta_diversity <- read.csv(file = "pairwase_beta.csv", header = T, row.names = 1)

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

p1 <- ggplot(data = epi_comparison, aes(x = Var1, y = Var2, fill = Beta_Diversity)) + geom_tile() + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred", midpoint = 0.13, name = "Weighted UniFrac") + labs(x = NULL, y = NULL, title ="Epilimnia") 


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

p2 <- ggplot(data = hypo_comparison, aes(x = Var1, y = Var2, fill = Beta_Diversity)) + geom_tile() + scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred", midpoint = 0.16, name = "Weighted UniFrac") + labs(x = NULL, y = NULL, title ="Hypolimnia")
plot_grid(p1, p2, nrow = 1)


#### p-values for richness tests

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

z <- pairwise.wilcox.test(epi.data$epi.obs, epi.data$epi.lakes, p.adjust.method = "bonferroni")
z <- melt(z$p.value)
write.csv(z, file = "epi_pvalues.csv")

z <- pairwise.wilcox.test(hypo.data$hypo.obs, hypo.data$hypo.lakes, p.adjust.method = "bonferroni")
z <- melt(z$p.value)
write.csv(z, file = "hypo_pvalues.csv")

#### Sampling coverage

lakekey <- substr(colnames(otu_table), start = 1, stop = 2)
datekey <- extract_date(colnames(otu_table))
coverage <- data.frame(lakekey, datekey)

ggplot(data = coverage, aes(x = datekey, y = lakekey)) + geom_point() + labs(x = NULL, y = "Lake", title = "Sampling Coverage")

#### Comparison of beta diversity metrics

# Use a reduced dataset - Trout Bog Hypo 2007
TBH07 <- prune_samples(sampledata$Layer == "H" & sampledata$Bog == "TB" & sampledata$Year == "07", alldata)
x <- UniFrac(TBH07, weighted = T, normalize = T)
beta_table <- melt(as.matrix(x))
colnames(beta_table) <- c("Sample1", "Sample2", "Weighted_UniFrac")

x <- UniFrac(TBH07, weighted = F, normalize = T)
unweighted <- melt(as.matrix(x))
beta_table$Unweighted_UniFrac <- unweighted$value

TBH07 <- bog_subset("TBH", otu_table)
TBH07 <- year_subset("07", TBH07)


x <- vegdist(t(TBH07), method = "bray")
bray <- melt(as.matrix(x))
beta_table$Bray_Curtis <- bray$value

x <- vegdist(t(TBH07), method = "jaccard")
jaccard <- melt(as.matrix(x))
beta_table$Jaccard <- jaccard$value

p1 <- ggplot(data = beta_table, aes(x = Weighted_UniFrac, y = Unweighted_UniFrac)) + geom_point() + labs(title = c(paste("r^2 =", round(cor(beta_table$Weighted_UniFrac, beta_table$Unweighted_UniFrac), 2))))
p2 <- ggplot(data = beta_table, aes(x = Weighted_UniFrac, y = Bray_Curtis)) + geom_point() + labs(title = c(paste("r^2 =", round(cor(beta_table$Weighted_UniFrac, beta_table$Bray_Curtis), 2))))
p3 <- ggplot(data = beta_table, aes(x = Weighted_UniFrac, y = Jaccard)) + geom_point() + labs(title = c(paste("r^2 =", round(cor(beta_table$Weighted_UniFrac, beta_table$Jaccard), 2))))
p4 <- ggplot(data = beta_table, aes(x = Unweighted_UniFrac, y = Bray_Curtis)) + geom_point() + labs(title = c(paste("r^2 =", round(cor(beta_table$Unweighted_UniFrac, beta_table$Bray_Curtis), 2))))
p5 <- ggplot(data = beta_table, aes(x = Unweighted_UniFrac, y = Jaccard)) + geom_point() + labs(title = c(paste("r^2 =", round(cor(beta_table$Unweighted_UniFrac, beta_table$Jaccard), 2))))
p6 <- ggplot(data = beta_table, aes(x = Bray_Curtis, y = Jaccard)) + geom_point() + labs(title = c(paste("r^2 =", round(cor(beta_table$Bray_Curtis, beta_table$Jaccard), 2))))

beta_comparison <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2)

# How does each metric perform in PCoA?
TBH07 <- prune_samples(sampledata$Layer == "H" & sampledata$Bog == "TB" & sampledata$Year == "07", alldata)
x <- UniFrac(TBH07, weighted = T, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 1, stop = 3))
wu_eigen <- as.data.frame(pcoa$eig)

x <- UniFrac(TBH07, weighted = F, normalize = T)
pcoa <- betadisper(x, substr(labels(x), start = 1, stop = 3))
uu_eigen <- as.data.frame(pcoa$eig)

TBH07 <- bog_subset("TBH", otu_table)
TBH07 <- year_subset("07", TBH07)

x <- vegdist(t(TBH07), method = "bray")
pcoa <- betadisper(x, substr(labels(x), start = 1, stop = 3))
bc_eigen <- as.data.frame(pcoa$eig)

x <- vegdist(t(TBH07), method = "jaccard")
pcoa <- betadisper(x, substr(labels(x), start = 1, stop = 3))
j_eigen <- as.data.frame(pcoa$eig)

wu_eigen <- wu_eigen/sum(wu_eigen) * 100
uu_eigen <- uu_eigen/sum(uu_eigen) * 100
bc_eigen <- bc_eigen/sum(bc_eigen) * 100
j_eigen <- j_eigen/sum(j_eigen) * 100

all_eigen <- c(wu_eigen[1:5, ], uu_eigen[1:5, ], bc_eigen[1:5, ], j_eigen[1:5, ])
eigen_names <- c(rep("Weighted_UniFrac", 5), rep("Unweighted_UniFrac", 5), rep("BrayCurtis", 5), rep("Jaccard", 5))
axis <- c(rep(c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5"), 4))

eigenvalues <- data.frame(eigen_names, axis, all_eigen)
colnames(eigenvalues) <- c("Metric", "Axis", "Variance")

ggplot(data = eigenvalues, aes(x = axis, y = Variance, fill = Metric)) + geom_bar(stat = "identity", position = "dodge") + labs(y = "% Variance Explained")

#### Beta diversity time decay plots

# Use matrix from task 1

site <- c("CBE", "FBE", "WSE", "NSE", "TBE", "SSE", "HKE", "MAE", "CBH", "FBH", "WSH", "NSH", "TBH", "SSH", "HKH", "MAH")

# For each site, make plots of beta diversity of tn:t0 and tn:tn-1

for(i in 1:length(site)){
  rows <- grep(site[i], rownames(beta_diversity))
  cols <- grep(site[i], colnames(beta_diversity))
  lake <- beta_diversity[rows,cols]
  datekey <- extract_date(colnames(lake))
  
  #Compare to time 0
  start <- lake[which(datekey == min(datekey)), ]
  if(dim(start)[1] > 1){
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

#### Measure variation at seasonal and site level

# Use weighted UniFrac matrix to measure mean pairwise distance within sites and years
# Make scatterplots with mean as line - geom_dotplot

beta_diversity <- read.csv(file = "pairwase_beta.csv", header = T, row.names = 1)
site <- c("CBE", "FBE", "WSE", "NSE", "TBE", "SSE", "HKE", "MAE", "CBH", "FBH", "WSH", "NSH", "TBH", "SSH", "HKH", "MAH")
melted_diversity <- melt(as.matrix(beta_diversity))
melted_diversity$Site1 <- substr(melted_diversity$Var1, start = 1, stop = 3)
melted_diversity$Site2 <- substr(melted_diversity$Var2, start = 1, stop = 3)
within_site <- melted_diversity[which(melted_diversity$Site1 == melted_diversity$Site2), ]

#Remove anything with the letter U
within_site <- within_site[grep("U", within_site$Site1, invert = T), ]
within_site$Site1 <- factor(within_site$Site1, levels = rev(c("CBE","CBH", "FBE", "FBH", "WSE", "WSH", "NSE", "NSH", "TBE", "TBH", "SSE", "SSH", "HKE", "HKH", "MAE", "MAH")))
pal <- rev(c("#a6cee3", "#a6cee3", "#1f78b4", "#1f78b4", "#b2df8a", "#b2df8a", "#33a02c", "#33a02c", "#fb9a99", "#fb9a99", "#e31a1c", "#e31a1c", "#fdbf6f", "#fdbf6f", "#ff7f00", "#ff7f00"))
ggplot(data = within_site, aes(x = Site1, y = value, fill = Site1)) + geom_violin() + coord_flip() + theme(legend.position = "none") + scale_fill_manual(values = pal) + labs(y = "Pairwise Weighted UniFrac Distance", x = NULL)
#report mean for each group
means <- c()
for(i in 1:length(site)){
  values <- within_site$value[which(within_site$Site1 == site[i])]
  means[i] <- mean(values)
}

#Test pairwise significance
pairwise.wilcox.test(within_site$value, within_site$Site1, p.adjust.method = "bonferroni")

#Repeat with only 2007 data
within_site$Year <- substr(within_site$Var1, start = 9, stop = 10)
vary07 <- within_site[which(within_site$Year == "07"), ]
ggplot(data = vary07, aes(x = Site1, y = value)) + geom_boxplot()
ggplot(data = vary07, aes(x = Site1, y = value)) + geom_violin()
pairwise.wilcox.test(vary07$value, vary07$Site1, p.adjust.method = "bonferroni")


# For every date in every lake, get the average temperature in both the epilimnion at hypolimnion

lake <- c("CBE", "FBE", "WSE", "NSE", "TBE", "SSE", "HKE", "MAE")
metalimnion <- c(1, 1, 2, 2, 2, 2, 2, 2)
depth <- c(2, 2, 4, 4, 6, 8, 18, 20)

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

mean(mean_meta2$Temperature[which(mean_meta2$Group.1 == "Hypo" & mean_meta2$Group.2 == "TB" & substr(as.character(mean_meta2$Group.3), start = 1, stop = 4) == "2009")])

# Add heatmap of temp to background of richness over time plot
meta2 <- metadata[,c(1,2,3,4)]
meta2$Layer <- substr(meta2$Sample_Name, start = 3, stop = 3)
meta2 <- meta2[which(meta2$Layer == "E"), ]
meta2$Site <- substr(meta2$Sample_Name, start = 1, stop = 2)
meta2$Date <- extract_date(meta2$Sample_Name)
meta2$Year <- substr(as.character(meta2$Date), start = 1, stop = 4)

ggplot(meta2[which(meta2$Site == "TB" & meta2$Year == "2007"), ], aes(x = Date, y = Depth, z = Temperature)) + geom_contour() + stat_contour(geom = "polygon", aes(fill = ..level..)) 
+ scale_y_reverse()

#Now add richness plot
###########
#Figure S2 - richness over time

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

TB_meta <- meta2[which(meta2$Site == "TB" & meta2$Year == "2007"), ]

figS2a <- 

ggplot(meta2[which(meta2$Site == "TB" & meta2$Year == "2007"), ], aes(x = Date, y = (-Depth + 6) * 54.667, z = Temperature)) + geom_contour() + stat_contour(geom = "polygon", aes(fill = ..level..)) + scale_fill_gradientn(colors = c("#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#fee08b", "#fdae61", "#f46d43", "#d53e4f")) + theme(panel.background = element_rect(fill = '#2166ac'))
ggplot() + geom_line(data = TB_richness, aes(x = date, y = richness), size = 1) + labs(title = "Trout Bog", x = NULL, y = "Observed Richness")

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

figS2b <- ggplot() + geom_line(data = NS_richness, aes(x = date, y = richness), size = 1.2) + labs(title = "North Sparkling Bog", x = NULL, y = "Observed Richness")  + geom_point(data = NS_richness[match(NSHmixes, NS_richness$date), ], aes(x = date, y = richness), size = 2, colour = "red")
figS2 <- plot_grid(figS2a, figS2b, align = "h", nrow = 2, labels = c("A", "B"))
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots_2017/FigS2.pdf", figS2, base_aspect_ratio = 1.5, base_height = 5)


#How does rarefaction level change richness?

otu_table <- read.csv("C:/Users/Alex/Documents/deblurred_nonrarified_otutable.csv", header = T, row.names = 1)
#data(taxonomy)
taxonomy <- read.csv("C:/Users/Alex/Documents/nonrarified_clean_taxonomy.csv", header = T, row.names = 1, colClasses = c("character"))
seqs2 <- read.dna("C:/Users/Alex/Documents/bog_repseqs_nonrarified_aligned.good.filter.fasta", format = "fasta", as.matrix = T)
otu_table <- otu_table[match(dimnames(seqs2)[[1]], rownames(otu_table)), ]
taxonomy <- taxonomy[match(rownames(otu_table), rownames(taxonomy)), ]


# Compare histograms of TB vs MA using data from nonrarified table, 5000 rarified, and 2500 rarified.

#Non-rarified (use proportions)
prop_table <- sweep(otu_table, 2, colSums(otu_table), '/')
MAH <- bog_subset("MAH", prop_table)
TBH <- bog_subset("TBH", prop_table)
CBH <- bog_subset("CBH", prop_table)
MAH.rich <- apply(MAH, 2, obs_richness)
TBH.rich <- apply(TBH, 2, obs_richness)
CBH.rich <- apply(CBH, 2, obs_richness)
lakekey <- c(rep("MAH", length(MAH.rich)), rep("TBH", length(TBH.rich)), rep("CBH", length(CBH.rich)))
richness <- c(MAH.rich, TBH.rich, CBH.rich)
plot.rich <- data.frame(lakekey, richness)

p1 <- ggplot(data = plot.rich, aes(x = richness, fill = lakekey)) + geom_density(alpha = 0.5) + labs(title = "Non-rarefied", x = "Observed Richness", y = "Density") + theme(legend.position = "none")

#Rarefy by 5000
rarefied1 <- otu_table[, which(colSums(otu_table) >= 5000)]
for(i in 1:dim(rarefied1)[2]){
  weighted_otus <- rep(rownames(rarefied1), rarefied1[, i])
  keep <- sample(weighted_otus, 5000, replace = F)
  newcounts <- table(keep)
  hits <- match(rownames(rarefied1), names(newcounts))
  newcolumn <- newcounts[hits]
  newcolumn[which(is.na(newcolumn) == T)] <- 0
  rarefied1[,i] <- newcolumn
}

MAH <- bog_subset("MAH", rarefied1)
TBH <- bog_subset("TBH", rarefied1)
CBH <- bog_subset("CBH", rarefied1)
MAH.rich <- apply(MAH, 2, obs_richness)
TBH.rich <- apply(TBH, 2, obs_richness)
CBH.rich <- apply(CBH, 2, obs_richness)
lakekey <- c(rep("MAH", length(MAH.rich)), rep("TBH", length(TBH.rich)), rep("CBH", length(CBH.rich)))
richness <- c(MAH.rich, TBH.rich, CBH.rich)
plot.rich <- data.frame(lakekey, richness)

p2 <- ggplot(data = plot.rich, aes(x = richness, fill = lakekey)) + geom_density(alpha = 0.5) + labs(title = "Rarefied to 5000", x = "Observed Richness", y = "Density") + theme(legend.position = "none")

#Rarefy by 2500
rarefied1 <- otu_table[, which(colSums(otu_table) >= 2500)]
for(i in 1:dim(rarefied1)[2]){
  weighted_otus <- rep(rownames(rarefied1), rarefied1[, i])
  keep <- sample(weighted_otus, 2500, replace = F)
  newcounts <- table(keep)
  hits <- match(rownames(rarefied1), names(newcounts))
  newcolumn <- newcounts[hits]
  newcolumn[which(is.na(newcolumn) == T)] <- 0
  rarefied1[,i] <- newcolumn
}

MAH <- bog_subset("MAH", rarefied1)
TBH <- bog_subset("TBH", rarefied1)
CBH <- bog_subset("CBH", rarefied1)
MAH.rich <- apply(MAH, 2, obs_richness)
TBH.rich <- apply(TBH, 2, obs_richness)
CBH.rich <- apply(CBH, 2, obs_richness)
lakekey <- c(rep("MAH", length(MAH.rich)), rep("TBH", length(TBH.rich)), rep("CBH", length(CBH.rich)))
richness <- c(MAH.rich, TBH.rich, CBH.rich)
plot.rich <- data.frame(lakekey, richness)

p3 <- ggplot(data = plot.rich, aes(x = richness, fill = lakekey)) + geom_density(alpha = 0.5) + labs(title = "Rarefied to 2500", x = "Observed Richness", y = "Density") + theme(legend.position = c(0.8, 0.8), legend.title = element_blank())


rarefaction <- plot_grid(p1, p2, p3, nrow = 1)


#Make Fig2 with different colorations. Add PERMANOVA table

# # Set up phyloseq object to run UniFrac on
data(metadata)
data(otu_table)
data(taxonomy)
seqs <- read.dna("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_repseqs_07Jul15.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(taxonomy))
sampledata <- sample_data(data.frame(Bog = substr(colnames(otu_table), start = 1, stop = 2), Layer = substr(colnames(otu_table), start = 3, stop = 3), Year = substr(colnames(otu_table), start = 9, stop = 10), row.names = colnames(otu_table), stringsAsfactors = F)) 
alldata <- phyloseq(OTU, TAX, sampledata, bogtree)

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
etalimnion <- c(1, 1, 2, 2, 2, 2, 2, 2)
depth <- c(2, 2, 4, 4, 6, 8, 18, 20)

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

ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_blank())
ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = as.integer(Date))) + geom_point(size = 1, alpha = 1/2) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + theme(legend.title = element_blank()) + scale_color_gradientn(colors = c("#ffffcc", "#a1dab4","#41b6c4", "#225ea8"))
ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Regime)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_blank())
ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Temp)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Epilimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_gradientn(colors = c("#2c7bb6", "#abd9e9","#fdae61", "#d7191c")) + theme(legend.title = element_blank())


adonis(x ~ lakes, as(sample_data(epi), "data.frame"))
# r^2 = 0.34498, p = 0.001
# Calculate PERMANOVA by mixing regime
adonis(x ~ regime, as(cbind(sample_data(epi), regime), "data.frame"))
# r^2 = 0.19574, p = 0.001

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

ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Lake)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_blank())
ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = as.integer(Date))) + geom_point(size = 1, alpha = 1/2) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + theme(legend.title = element_blank()) + scale_color_gradientn(colors = c("#ffffcc", "#a1dab4","#41b6c4", "#225ea8"))
ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Regime)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_blank())
ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Temp)) + geom_point(size = 1, alpha = 1/2) + labs(title = "Hypolimnia", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")")) + scale_color_gradientn(colors = c("#2c7bb6", "#abd9e9","#fdae61", "#d7191c")) + theme(legend.title = element_blank())



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
save_plot("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots_2017/Fig2.pdf", fig2, base_aspect_ratio = 1.2, base_height = 6)
