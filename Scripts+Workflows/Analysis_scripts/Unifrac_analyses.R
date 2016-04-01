library(OTUtable)
library(phyloseq)
library(ggplot2)

data(otu_table)
data(taxonomy)

# Unfortunately, the tree of all OTUs is too large for ape to write to file (recursion errors)
# A new tree has to be made each time the script is run. This takes approx. 15 min for our 6,208 sequences. I'd recommend making yourself a cup of tea during this time.
seqs <- read.dna("C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_repseqs_07Jul15.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

# Make OTU table, taxonomy, and sampledata datasets
OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(taxonomy))

sampledata <- sample_data(data.frame(Bog = substr(colnames(otu_table), start = 1, stop = 2), Layer = substr(colnames(otu_table), start = 3, stop = 3), Year = substr(colnames(otu_table), start = 9, stop = 10), row.names = colnames(otu_table), stringsAsfactors = F))                                                                                                                                          

alldata <- phyloseq(OTU, TAX, sampledata, bogtree)

TBH <- prune_samples(sampledata$Bog == "TB" & sampledata$Layer == "H", alldata)
TBH_year <- factor(substr(sample_names(TBH), start = 9, stop = 10), levels = c("05", "07", "08", "09"))

x <- UniFrac(TBH, weighted = T, normalize = T)
pcoa <- betadisper(x, TBH_year)
scores <- scores(pcoa)
TBHcentroids <- scores$centroids

plot.pcoa <- data.frame(scores$sites, TBH_year)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

p1 <- ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point(size=3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 1, size = 12, colour = "black"), axis.title = element_text(size = 15, vjust=0.7), axis.text.y = element_text(colour="black", size=12))

adonis(x ~ Year, as(sample_data(TBH_physeq), "data.frame"))
