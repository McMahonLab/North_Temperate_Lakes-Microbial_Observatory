#Convert OTU table to count table containing only a subset of samples. Adjust fasta file to match sequences that still contain reads.

#Load table
library(OTUtable)
data(otu_table)

args <- commandArgs(trailingOnly=TRUE)
#Subset samples wanted
bog <- args[1]

bog_table <- bog_subset(bog, otu_table)

#Add rep seqs and totals column
repseqs <- rownames(bog_table)
repseqs2 <- paste("bog", repseqs, collapse=NULL, sep="")
totals <- rowSums(bog_table)

countfile <- data.frame(Representative_Sequences=repseqs2, totals=totals, bog_table)
countfile <- countfile[which(countfile$totals > 0),]
#Output file. Turn off rownames so that we don't get duplicate columns
write.table(countfile, file=paste("C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/UniFrac_analysis/Annual_variation/temp.count", sep = ""), quote=F, sep="\t", row.names=F, col.names=T)

#Keep only sequences in the rep seq fasta file that are in the new count file
fasta <- read.table("C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/UniFrac_analysis/bog_repseqs_22Jan16_10charnames.fasta", header=F, colClasses = c("character"))

#sequences to keep
otuIDs <- as.character(countfile$Representative_Sequences)
#Spit IDs and sequences into separate vectors
fastaIDs <- fasta$V1[seq(1, by=2, len=length(fasta$V1))]
fastaseqs <- fasta$V1[seq(2, by=2, len=length(fasta$V1))]
fastaIDs <- fastaIDs[which(is.na(fastaIDs) == F)]
fastaseqs <- fastaseqs[which(is.na(fastaseqs) == F)]
#Remove ">" from fastaIDs
fastaIDs <- substr(fastaIDs, start=2, stop=11)
#Match sequences to keep and fastaIDs
keep <- match(otuIDs, fastaIDs)
keep.fastaIDs <- fastaIDs[keep]
keep.fastaseqs <- fastaseqs[keep]
#Add carrot back onto the fasta IDS
keep.fastaIDs <- paste(">", keep.fastaIDs, sep = "")
#Interleave the two fasta vectors
new.fasta <- c(rbind(keep.fastaIDs, keep.fastaseqs))

write.table(new.fasta, file="C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/UniFrac_analysis/Annual_variation/temp.fasta", quote=F, row.names=F, col.names=F)
