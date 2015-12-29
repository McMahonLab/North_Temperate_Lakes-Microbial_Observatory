library(OTUtable)
library(indicspecies)

data(otu_table)
data(taxonomy)

tribe_table <- combine_otus("Tribe", otu_table, taxonomy)
clade_table <- combine_otus("Clade", otu_table, taxonomy)
lineage_table <- combine_otus("Lineage", otu_table, taxonomy)
order_table <- combine_otus("Order", otu_table, taxonomy)
class_table <- combine_otus("Class", otu_table, taxonomy)
phylum_table <- combine_otus("Phylum", otu_table, taxonomy)

named_otu_table <- otu_table
fullnames <- c()
for(i in 1:dim(taxonomy)[1]){
  fullnames[i] <- paste(taxonomy[i,], collapse = ";")
}
fullnames <- make.unique(fullnames)
rownames(named_otu_table) <- fullnames

full_table <- rbind(named_otu_table, tribe_table, clade_table, lineage_table, order_table, class_table, phylum_table)

classified <- grep("unclassified", rownames(full_table))
classified1 <- grep("__$", rownames(full_table))
input_table <- full_table[-c(classified, classified1),]

threshold <- quantile(rowSums(input_table))[4]
input_table <- input_table[which(rowSums(input_table) >= threshold),]
input_table <- t(input_table)
input_table <- as.data.frame(input_table)

lakeid <- c("CB", "FB", "WS", "NS", "TB", "SS", "HK","MA")

lakes <- substr(rownames(input_table), start=1, stop=2)

mixing_groups <- c()
mixing_groups[which(lakes == "CB" | lakes == "FB" | lakes == "WS")] <- 1
mixing_groups[which(lakes == "TB"| lakes == "NS"| lakes == "SS")] <- 2
mixing_groups[which(lakes == "MA" | lakes == "HK")] <- 3

clades_by_mixing <- multipatt(x = input_table, cluster = mixing_groups, func = "r.g", control = how(nperm = 9999))

write.csv(clades_by_mixing$sign, file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/indicators_by_mixing.csv")
