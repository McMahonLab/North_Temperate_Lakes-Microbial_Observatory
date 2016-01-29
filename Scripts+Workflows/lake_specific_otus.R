#how many OTUs are lake specific?

library(OTUtable)
data(otu_table)
data(taxonomy)

#Look just at Actinobacterial OTUs
actinos <- grab_group("acV-A1", "Tribe", otu_table, taxonomy)
actinos <- actinos[which(rowSums(actinos) > 500),]
CBH <- bog_subset("CBH", actinos)
TBH <- bog_subset("TBH", actinos)
NSH <- bog_subset("NSH", actinos)
SSH <- bog_subset("SSH", actinos)
MAH <- bog_subset("MAH", actinos)

#Make heatmap of acI spread

totals <- cbind(rowSums(CBH), rowSums(NSH), rowSums(TBH), rowSums(SSH), rowSums(MAH))
#Remove bacI-B
totals <- totals[c(2:8),]
colnames(totals) <- c("CBH", "NSH", "TBH", "SSH", "MAH")
totals <- reduce_names(totals)
heatmap(totals, Rowv=NA, Colv=NA)
