library(OTUtable)
library(reshape2)
library(dplyr)
library(ggplot2)
library(grid)
library(raster)
library(ggrepel)
library(ape)
library(phyloseq)

data(metadata)

seq_table <- read.csv("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_seqstable.csv", row.names = 1)
seq_taxonomy <- read.csv(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/seqs.98.cleantaxonomy.csv", row.names = 1, header = T, colClasses = c("character"))

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



clade_table <- combine_otus("Clade", seq_table, seq_taxonomy)

clades <- clade_table[grep("g__|unclassified", rownames(clade_table), invert = T), ]
#Keep only the clade name
# I'm only showing a few lakes, but these can be adjusted easilty
lakes <- c("CBE", "TBE", "SSE", "MAE")

# For each lake, calculate my 3 metrics and plot. Save the plot to combine  into 1 pdf
for(i in 1:length(lakes)){
  table <- bog_subset(lakes[i], clades)
  #Remove clades with very few reads
  if(i == 4){
    table <- table[which(rowSums(table) > 400),]
  }else{
    table <- table[which(rowSums(table) > 10),]
  }
  
  
  abundance <- c()
  persistance <- c()
  variance <- c()
  
  for(j in 1:dim(table)[1]){
    row <- table[j, ]
    abundance[j] <- sum(row[which(row > 0)])/length(which(row > 0))
    persistence[j] <- length(which(row > 0))/length(row)
    variance[j] <- cv(as.numeric(row))
  }
  
  to.plot <- data.frame(abundance, persistence, variance)
  to.plot$clade <- rownames(reduce_names(table))
  
  z <- ggplot(data = to.plot, aes(x = variance, y = abundance/2500, fill = persistence)) + geom_point() + theme_bw() + labs(title = paste(lakes[i]), x = "Variability (CV)", y = "Mean Abundance when Present") + geom_label_repel(aes(label = clade), size = 2, force = 12, box.padding = unit(0.01, "lines")) + scale_fill_gradient(low = "white", high = "lightgreen")
  
  assign(paste("p", i, sep = ""), z)
}

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figure5.pdf", width = 3.3125*3, height = 9)
multiplot(p1, p2, p3, p4, cols = 2)
dev.off()

year <- c("05", "07", "08", "09")
for(i in 1:length(year)){
  table <- bog_subset("TBH", clades)
  table <- year_subset(year[i], table)
  table <- table[which(rowSums(table) > 10),]
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
  to.plot$clade <- rownames(reduce_names(table))
  
  z <- ggplot(data = to.plot, aes(x = variance, y = abundance/2500, fill = persistance)) + geom_point() + theme_bw() + labs(title = paste(year[i]), x = "Variability (CV)", y = "Mean Abundance when Present") + geom_label_repel(aes(label = clade), size = 2, force = 12, box.padding = unit(0.01, "lines")) + scale_fill_gradient(low = "white", high = "lightgreen")
  
  assign(paste("p", i, sep = ""), z)
}

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS7.pdf", width = 3.3125*3, height = 9)
multiplot(p1, p2, p3, p4, cols = 2)
dev.off()

seqs <- read.dna("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/qc.bogs.clean.min25.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

# Set up phyloseq object to run UniFrac on
OTU <- otu_table(as.matrix(seq_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(seq_taxonomy))
sampledata <- sample_data(data.frame(Bog = substr(colnames(seq_table), start = 1, stop = 2), Layer = substr(colnames(seq_table), start = 3, stop = 3), Year = substr(colnames(seq_table), start = 9, stop = 10), row.names = colnames(seq_table), stringsAsfactors = F)) 
alldata <- phyloseq(OTU, TAX, sampledata, bogtree)

# Select lake and layer, then calculate UniFrac and plot vs time difference
TBH <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "H", alldata)
x <- UniFrac(TBH, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
colors[which(x$Mon1 == "NOV" & x$Mon2 == "NOV")] <- "blue"
p1 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Trout Bog Hypolimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

TBE <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "E", alldata)
x <- UniFrac(TBE, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
x$Mon1 <- substr(x$Var1, start = 6, stop = 8)
x$Mon2 <- substr(x$Var2, start = 6, stop = 8)
colors <- rep("black", dim(x)[1])
p2 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Trout Bog Epilimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

MAH <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer == "H", alldata)
x <- UniFrac(MAH, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
x <- x[which(x$Var1 != "MAH23OCT07.R1" & x$Var1 != "MAH30JUL07.R1" & x$Var2 != "MAH23OCT07.R1" & x$Var2 != "MAH30JUL07.R1" & x$Var1 != "MAH23OCT07.R2" & x$Var1 != "MAH30JUL07.R2" & x$Var2 != "MAH23OCT07.R2" & x$Var2 != "MAH30JUL07.R2"),]
p3 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Mary Lake Hypolimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figure3.pdf", width = 3.3125*2, height = 7)
multiplot(p1, p2, p3, cols = 1)
dev.off()

# Select lake and layer, then calculate UniFrac and plot vs time difference
NSH <- prune_samples(sampledata$Bog == "NS" & sampledata$Layer == "H", alldata)
x <- UniFrac(NSH, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p1 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "North Sparkling Hypolimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

NSE <- prune_samples(sampledata$Bog == "NS" & sampledata$Layer == "E", alldata)
x <- UniFrac(NSE, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p2 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "North Sparkling Epilimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

MAE <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer == "E", alldata)
x <- UniFrac(MAE, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p3 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Mary Lake Epilimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS8.1.pdf", width = 3.3125*2, height = 7)
multiplot(p1, p2, p3, cols = 1)
dev.off()

SSH <- prune_samples(sampledata$Bog == "SS" & sampledata$Layer == "H", alldata)
x <- UniFrac(SSH, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p1 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "South Sparkling Hypolimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

SSE <- prune_samples(sampledata$Bog == "SS" & sampledata$Layer == "E", alldata)
x <- UniFrac(SSE, weighted = T, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
p2 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "South Sparkling Epilimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()


pdf(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Plots/figureS8.2.pdf", width = 3.3125*2, height = 7/3*2)
multiplot(p1, p2, cols = 1)
dev.off()

#Futz with the UniFrac parameters
TBH <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "H", alldata)
TBH <- prune_taxa(taxa_sums(TBH) > 1000, TBH)
x <- UniFrac(TBH, weighted = F, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
colors[which(x$Mon1 == "NOV" & x$Mon2 == "NOV")] <- "blue"
p1 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Trout Bog Hypolimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

TBE <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "E", alldata)
TBE <- prune_taxa(taxa_sums(TBE) > 2000, TBE)
x <- UniFrac(TBE, weighted = F, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
x$Mon1 <- substr(x$Var1, start = 6, stop = 8)
x$Mon2 <- substr(x$Var2, start = 6, stop = 8)
colors <- rep("black", dim(x)[1])
#p2 <- 
  ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Trout Bog Epilimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

MAH <- prune_samples(sampledata$Bog == "MA" & sampledata$Layer == "H", alldata)
MAH <- prune_taxa(taxa_sums(MAH) > 1000, MAH)
x <- UniFrac(MAH, weighted = F, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
x <- x[which(x$Var1 != "MAH23OCT07.R1" & x$Var1 != "MAH30JUL07.R1" & x$Var2 != "MAH23OCT07.R1" & x$Var2 != "MAH30JUL07.R1" & x$Var1 != "MAH23OCT07.R2" & x$Var1 != "MAH30JUL07.R2" & x$Var2 != "MAH23OCT07.R2" & x$Var2 != "MAH30JUL07.R2"),]
p3 <- ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Mary Lake Hypolimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

multiplot(p1, p2, p3, cols = 1)


#Compare UniFrac to Bray-Curtis

TBE <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "E", alldata)
TBE <- prune_taxa(taxa_sums(TBE) > 2000, TBE)
TBE2 <- TBE@otu_table@.Data
TBE3 <- combine_otus("Phylum", TBE2, seq_taxonomy[match(rownames(TBE2), rownames(seq_taxonomy)),])
x <- vegdist(t(TBE3), method = "bray")
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Trout Bog Epilimnion", x = "Time Difference", y = "Bray Curtis") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()


# Is there a sweet spot of correlation vs genetic distance?
# Reduce dataset to actinos only in TBE

TBE_actinos <- bog_subset("TBE", seq_table)
TBE_actinos <- TBE_actinos[match(rownames(seq_taxonomy)[which(seq_taxonomy$Phylum == "p__Actinobacteria")], rownames(TBE_actinos)),]
ac_index <- match(rownames(seq_taxonomy)[which(seq_taxonomy$Phylum == "p__Actinobacteria")], rownames(d))
dist_actinos <- d[ac_index, ac_index]

#Correlate all of the actinos to each other
y <- cor(t(TBE_actinos))
y <- melt(y)
z <- melt(dist_actinos)
y$dist <- z$value
y <- y[which(is.na(y$value) == F),]
y <- y[which(y$value != 1),]
y$dist <- as.factor(y$dist)


#That's not helpful. No change that I can see

#Focus now on fitting a sine wave or logarithmic model to the regression

TBE <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "E", alldata)
x <- UniFrac(TBE, weighted = F, normalize = T)
x <- melt(as.matrix(x))
x$Date <- as.numeric(abs(extract_date(x$Var1) - extract_date(x$Var2)))
x$Mon1 <- substr(x$Var1, start = 6, stop = 8)
x$Mon2 <- substr(x$Var2, start = 6, stop = 8)
colors <- rep("black", dim(x)[1])
#p2 <- 
ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Trout Bog Epilimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + geom_smooth()

#Fit the sine wave equation
wave <- lm(value ~ sin(Date) + cos(Date), x)
ggplot(x, aes(x = Date, y = value)) + geom_point(size = 0.5, alpha = 1/10) + theme_bw() + labs(title = "Trout Bog Epilimnion", x = "Time Difference", y = "UniFrac Distance") + scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4)) + stat_function(fun = function(Date) 0.25*(sin(2*pi*(Date+365)/365) + cos(2*pi*(Date)/365)) )

ssp <- spectrum(x$value)  
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
reslm <- lm(x$value ~ sin(2*pi/per*x$Date)+cos(2*pi/per*x$Date))
summary(reslm)

rg <- diff(range(x$value))
plot(x$value~x$Date,ylim=c(min(x$value)-0.1*rg,max(x$value)+0.1*rg))
lines(fitted(reslm)~x$Date,col=4,lty=2)   # dashed blue line is sin fit

# including 2nd harmonic really improves the fit
reslm2 <- lm(x$value ~ sin(2*pi/365*x$Date)+cos(2*pi/365*x$Date)+sin(4*pi/365*x$Date)+cos(4*pi/365*x$Date))
summary(reslm2)
lines(fitted(reslm2)~x$Date,col=3)
