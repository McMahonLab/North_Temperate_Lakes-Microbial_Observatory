# Compare unweighted + weighted UniFrac vs Bray Curtis

# Step 1: Type II error (type I already tested and is negligable for all 3 metrics)

library(OTUtable)
library(vegan)
library(phyloseq)
library(ape)

data(otu_table)
data(taxonomy)
data(metadata)

seqs <- read.dna("C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_repseqs_07Jul15.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

MAHtable <- bog_subset("MAH", otu_table)
MAHtable07 <- year_subset("07", MAHtable)
REP07 <- rowSums(MAHtable07)/dim(MAHtable07)[2]

test <- seq(from = 250, to = 999, by = 15)

fakesamples07 <- matrix(ncol = 50, nrow = 6208)
for(i in 1:length(test)){
  fakesample <- c(rep(0, 6208))
  matched.otus <- sample(which(REP07 > 0), size = test[i],  replace = F)
  fakesample[matched.otus] <- REP07[matched.otus]
  remaining.counts <- 2500 - sum(fakesample)
  other.otus <- sample(rownames(otu_table), size = remaining.counts, prob = rowSums(otu_table), replace = T)
  fake.counts <- table(other.otus)
  fakesample[match(names(fake.counts), names(REP07))] <- fakesample[match(names(fake.counts), names(REP07))] + fake.counts
  fakesample[which(is.na(fakesample) == T)] <- 0
  fakesamples07[,i] <- fakesample
}
fakesamples07 <- data.frame(fakesamples07, REP07)
rownames(fakesamples07) <- rownames(otu_table)

# Make phyloseq object
OTU <- otu_table(as.matrix(fakesamples07), taxa_are_rows = T)
TAX <- tax_table(as.matrix(taxonomy))
testdata07 <- phyloseq(OTU, TAX, bogtree)

x <- UniFrac(testdata07, weighted = T, normalize = T)
x <- as.matrix(x)
UFW07 <- x[1:50,51]

y <- UniFrac(testdata07, weighted = F, normalize = T)
y <- as.matrix(y)
UF07 <- y[1:50,51]

z <- vegdist(t(fakesamples07), method = "bray")
z <- as.matrix(z)
BC07 <- z[1:50,51]

UFW.typeII <- c()
UF.typeII <- c()
BC.typeII <- c()

for(k in 1:10000){
  shuffle <- sample(1:50, size = 50, replace = F)
  model <- lm(UFW07 ~ shuffle)
  UFW.typeII[k] <- summary(model)$coefficients[2,4]
  
  model <- lm(UF07 ~ shuffle)
  UF.typeII[k] <- summary(model)$coefficients[2,4]
  
  model <- lm(BC07 ~ shuffle)
  BC.typeII[k] <- summary(model)$coefficients[2,4]
}

length(which(UFW.typeII < 0.05))/1000
length(which(UF.typeII < 0.05))/1000
length(which(BC.typeII < 0.05))/1000


#Ran this multiple times - all seem roughly equilivalent in terms of type II error

#Try the PCoA variation test
#Do this both on TBH by year and TBH vs MAH 07

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

OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(taxonomy))

sampledata <- sample_data(data.frame(Bog = substr(colnames(otu_table), start = 1, stop = 2), Layer = substr(colnames(otu_table), start = 3, stop = 3), Year = substr(colnames(otu_table), start = 9, stop = 10), row.names = colnames(otu_table), stringsAsfactors = F))                                                                                                                                      
alldata <- phyloseq(OTU, TAX, sampledata, bogtree)
years <- c("05", "07", "08", "09")
colors <- c("#a6cee3", "#1f78b4", "#33a02c", "#b2df8a")

TBH <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "H", alldata)
TBH_year <- factor(substr(sample_names(TBH), start = 9, stop = 10), levels = years)

x1 <- UniFrac(TBH, weighted = F, normalize = T)
pcoa1 <- betadisper(x1, TBH_year)
scores1 <- scores(pcoa1)
plot.pcoa1 <- data.frame(scores1$sites, TBH_year)
colnames(plot.pcoa1) <- c("PCoA1", "PCoA2", "Year")
uw1 <- round(pcoa1$eig[1]/sum(pcoa1$eig), digits = 2)
uw2 <- round(pcoa1$eig[2]/sum(pcoa1$eig), digits = 2)

p1 <- ggplot(data=plot.pcoa1, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + labs(title = "unweighted UniFrac", x = paste("PCoA1 (", uw1, ")", sep = ""), y = paste("PCoA2 (", uw2, ")", sep = "")) + scale_color_manual(values = colors) 

x2 <- UniFrac(TBH, weighted = T, normalize = T)
pcoa2 <- betadisper(x2, TBH_year)
scores2 <- scores(pcoa2)
plot.pcoa2 <- data.frame(scores2$sites, TBH_year)
colnames(plot.pcoa2) <- c("PCoA1", "PCoA2", "Year")
w1 <- round(pcoa2$eig[1]/sum(pcoa2$eig), digits = 2)
w2 <- round(pcoa2$eig[2]/sum(pcoa2$eig), digits = 2)

p2 <- ggplot(data=plot.pcoa2, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + labs(title = "weighted UniFrac", x = paste("PCoA1 (", w1, ")", sep = ""), y = paste("PCoA2 (", w2, ")", sep = "")) + scale_color_manual(values = colors) 

x3 <- distance(TBH, method = "bray", type = "samples")
pcoa3 <- betadisper(x3, TBH_year)
scores3 <- scores(pcoa3)
plot.pcoa3 <- data.frame(scores3$sites, TBH_year)
colnames(plot.pcoa3) <- c("PCoA1", "PCoA2", "Year")
bc1 <- round(pcoa3$eig[1]/sum(pcoa3$eig), digits = 2)
bc2 <- round(pcoa3$eig[2]/sum(pcoa3$eig), digits = 2)

p3 <- ggplot(data=plot.pcoa3, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + labs(title = "Bray-Curtis", x = paste("PCoA1 (", bc1, ")", sep = ""), y = paste("PCoA2 (", bc2, ")", sep = "")) + scale_color_manual(values = colors) 

multiplot(p1, p2, p3, cols = 3)

x4 <- distance(TBH, method = "jaccard", type = "samples")
pcoa4 <- betadisper(x4, TBH_year)
scores4 <- scores(pcoa4)
plot.pcoa4 <- data.frame(scores4$sites, TBH_year)
colnames(plot.pcoa4) <- c("PCoA1", "PCoA2", "Year")
j1 <- round(pcoa4$eig[1]/sum(pcoa4$eig), digits = 2)
j2 <- round(pcoa4$eig[2]/sum(pcoa4$eig), digits = 2)

p4 <- ggplot(data=plot.pcoa4, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position="none") + labs(title = "Jaccard", x = paste("PCoA1 (", j1, ")", sep = ""), y = paste("PCoA2 (", j2, ")", sep = "")) + scale_color_manual(values = colors) 
#Looks like weighted UniFrac wins in this case

#Try the TBH vs MAH 07

#Make PCoA of TB vs MA
TBH_MAH <- prune_samples(sampledata$Year == "07" & sampledata$Bog == "TB" & sampledata$Layer == "H" | sampledata$Year == "07"  & sampledata$Bog == "MA" & sampledata$Layer == "H", alldata)
TBH_MAH_lake <- factor(substr(sample_names(TBH_MAH), start = 1, stop = 2), levels = c("TB", "MA"))
TBH_MAH_dates <- extract_date(sample_names(TBH_MAH))
even <- estimate_richness(TBH_MAH, measures = "Shannon")/log(length(which(taxa_sums(TBH_MAH) > 0)))

x1 <- UniFrac(TBH_MAH , weighted = F, normalize = T)
pcoa1 <- betadisper(x1, TBH_MAH_lake)

scores1 <- scores(pcoa1)

plot.pcoa1 <- data.frame(scores1$sites, TBH_MAH_lake, TBH_MAH_dates, even)
colnames(plot.pcoa1) <- c("PCoA1", "PCoA2", "Lake", "Date", "Evenness")
plot.pcoa1 <- plot.pcoa1[order(plot.pcoa1$Date),]

#Calculate weights of axes
uw1 <- round(pcoa1$eig[1]/sum(pcoa1$eig), digits = 2)
uw2 <- round(pcoa1$eig[2]/sum(pcoa1$eig), digits = 2)

p1 <- ggplot(data=plot.pcoa1, aes(x = PCoA1, y = PCoA2, color = Evenness, shape = Lake)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_color_gradientn(colors = c("firebrick3", "gold2", "springgreen4")) + geom_path(data = plot.pcoa1[which(plot.pcoa1$Lake == "TB"),], aes(x = PCoA1, y = PCoA2), size = 0.25) + geom_path(data = plot.pcoa1[which(plot.pcoa1$Lake == "MA"),], aes(x = PCoA1, y = PCoA2), size = 0.25) + labs(x = paste("PCoA1 (", uw1, ")", sep = ""),  y = paste("PCoA2 (", uw2, ")", sep = ""), title = "unweighted UniFrac")

x2 <- UniFrac(TBH_MAH , weighted = T, normalize = T)
pcoa2 <- betadisper(x2, TBH_MAH_lake)

scores2 <- scores(pcoa2)

plot.pcoa2 <- data.frame(scores2$sites, TBH_MAH_lake, TBH_MAH_dates, even)
colnames(plot.pcoa2) <- c("PCoA1", "PCoA2", "Lake", "Date", "Evenness")
plot.pcoa2 <- plot.pcoa2[order(plot.pcoa2$Date),]

#Calculate weights of axes
w1 <- round(pcoa2$eig[1]/sum(pcoa2$eig), digits = 2)
w2 <- round(pcoa2$eig[2]/sum(pcoa2$eig), digits = 2)

p2 <- ggplot(data=plot.pcoa2, aes(x = PCoA1, y = PCoA2, color = Evenness, shape = Lake)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_color_gradientn(colors = c("firebrick3", "gold2", "springgreen4")) + geom_path(data = plot.pcoa2[which(plot.pcoa2$Lake == "TB"),], aes(x = PCoA1, y = PCoA2), size = 0.25) + geom_path(data = plot.pcoa2[which(plot.pcoa2$Lake == "MA"),], aes(x = PCoA1, y = PCoA2), size = 0.25) + labs(x = paste("PCoA1 (", w1, ")", sep = ""),  y = paste("PCoA2 (", w2, ")", sep = ""), title = "weighted UniFrac")

x3 <- distance(TBH_MAH, method = "bray", type = "samples")
pcoa3 <- betadisper(x3, TBH_MAH_lake)

scores3 <- scores(pcoa3)

plot.pcoa3 <- data.frame(scores3$sites, TBH_MAH_lake, TBH_MAH_dates, even)
colnames(plot.pcoa3) <- c("PCoA1", "PCoA2", "Lake", "Date", "Evenness")
plot.pcoa3 <- plot.pcoa3[order(plot.pcoa3$Date),]

#Calculate weights of axes
bc1 <- round(pcoa3$eig[1]/sum(pcoa3$eig), digits = 2)
bc2 <- round(pcoa3$eig[2]/sum(pcoa3$eig), digits = 2)

p3 <- ggplot(data=plot.pcoa3, aes(x = PCoA1, y = PCoA2, color = Evenness, shape = Lake)) + geom_point(size=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 12, colour = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), axis.text.y = element_text(colour = "black", size = 10), panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_color_gradientn(colors = c("firebrick3", "gold2", "springgreen4")) + geom_path(data = plot.pcoa3[which(plot.pcoa3$Lake == "TB"),], aes(x = PCoA1, y = PCoA2), size = 0.25) + geom_path(data = plot.pcoa3[which(plot.pcoa3$Lake == "MA"),], aes(x = PCoA1, y = PCoA2), size = 0.25) + labs(x = paste("PCoA1 (", bc1, ")", sep = ""),  y = paste("PCoA2 (", bc2, ")", sep = ""), title = "Bray-Curtis")

multiplot(p1, p2, p3, cols = 3)
