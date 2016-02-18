#Set up environment
library(OTUtable)
library(indicspecies)
data(otu_table)
data(taxonomy)

#Make tables at each phylogenetic level
tribe_table <- combine_otus("Tribe", otu_table, taxonomy)
clade_table <- combine_otus("Clade", otu_table, taxonomy)
lineage_table <- combine_otus("Lineage", otu_table, taxonomy)
order_table <- combine_otus("Order", otu_table, taxonomy)
class_table <- combine_otus("Class", otu_table, taxonomy)
phylum_table <- combine_otus("Phylum", otu_table, taxonomy)

#Change OTU number designations to full taxonomic assignment
named_otu_table <- otu_table
fullnames <- c()
for(i in 1:dim(taxonomy)[1]){
  fullnames[i] <- paste(taxonomy[i,], collapse = ";")
}
fullnames <- make.unique(fullnames)
rownames(named_otu_table) <- fullnames

#Combine tables at each level into one giant table
full_table <- rbind(named_otu_table, tribe_table, clade_table, lineage_table, order_table, class_table, phylum_table)

#Remove groups that are unclassified at any level
classified <- grep("unclassified", rownames(full_table))
classified1 <- grep("__$", rownames(full_table))
input_table <- full_table[-c(classified, classified1),]

#Keep only the top quantile in abundance
threshold <- quantile(rowSums(input_table))[4]
input_table <- input_table[which(rowSums(input_table) >= threshold),]

#Format for input into mulitplatt() function
input_table <- t(input_table)
input_table <- as.data.frame(input_table)

#Group by mixing regime
lakeid <- c("CB", "FB", "WS", "NS", "TB", "SS", "HK","MA")
lakes <- substr(rownames(input_table), start=1, stop=2)

mixing_groups <- c()
mixing_groups[which(lakes == "CB" | lakes == "FB" | lakes == "WS")] <- 1
mixing_groups[which(lakes == "TB"| lakes == "NS"| lakes == "SS")] <- 2
mixing_groups[which(lakes == "MA" | lakes == "HK")] <- 3

#Run indicator taxa analysis
clades_by_mixing <- multipatt(x = input_table, cluster = mixing_groups, func = "r.g", control = how(nperm = 9999))

#Save data
#write.csv(clades_by_mixing$sign, file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/indicators_by_mixing.csv")

#Get top indictators of polymictic lakes (group 1)
results <- clades_by_mixing$sign
total <- sum(rowSums(full_table))
poly <- results[which(results$index == 1),]
poly <- poly[order(poly$stat, decreasing=T),]

poly.indicators <- poly[c(1, 3, 4, 5, 11, 12, 13, 18, 25, 26), 4:5]
find <- match(rownames(poly.indicators), rownames(full_table))
abundance <- rowSums(full_table[find,])/total
poly.indicators$Abundance <- abundance
poly.indicators$p.value <- NULL
poly.indicators$index <- NULL
colnames(poly.indicators) <- c("Phi", "% of Community")
write.csv(poly.indicators, file = "C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/polymictic_indicators.csv", row.names=T)

dim <- results[which(results$index == 2),]
dim <- dim[order(dim$stat, decreasing=T),]

dim.indicators <- dim[c(1, 6, 9, 11, 12, 13, 14, 15, 16, 20), 4:5]
find <- match(rownames(dim.indicators), rownames(full_table))
abundance <- rowSums(full_table[find,])/total
dim.indicators$Abundance <- abundance
dim.indicators$p.value <- NULL
dim.indicators$index <- NULL
colnames(dim.indicators) <- c("Phi", "% of Community")
write.csv(dim.indicators, file = "C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/dimictic_indicators.csv", row.names=T)

mero <- results[which(results$index == 3),]
mero <- mero[order(mero$stat, decreasing=T),]

mero.indicators <- mero[c(1, 7, 8, 10, 16, 17, 22, 25, 27, 29), 4:5]
find <- match(rownames(mero.indicators), rownames(full_table))
abundance <- rowSums(full_table[find,])/total
mero.indicators$Abundance <- abundance
mero.indicators$p.value <- NULL
mero.indicators$index <- NULL
colnames(mero.indicators) <- c("Phi", "% of Community")
write.csv(mero.indicators, file = "C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/meromictic_indicators.csv", row.names=T)
