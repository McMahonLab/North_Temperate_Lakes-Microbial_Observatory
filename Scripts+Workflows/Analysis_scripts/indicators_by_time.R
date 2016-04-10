library(OTUtable)
library(indicspecies)
data(otu_table)
data(taxonomy)

#Keep only TBH
TBH_table <- bog_subset("TBH", otu_table)

#Keep only OTUs that appear at least 10 times with at least 
persistent <- c()
abundant <- c()
for(i in 1:dim(TBH_table)[1]){
  group <- TBH_table[i,]
  persistent[i] <- length(which(group > 0)) >= 10
  abundant[i] <- sum(group) >= 25
}

TBH_table <- TBH_table[which(persistent == T & abundant == T), ]

TBH_taxonomy <- taxonomy[which(rownames(taxonomy) %in% rownames(TBH_table) == T),]

# Change OTU number designations to full taxonomic assignment
fullnames <- c()
for(i in 1:dim(TBH_taxonomy)[1]){
  fullnames[i] <- paste(TBH_taxonomy[i, ], collapse = ";")
}
fullnames <- make.unique(fullnames)
rownames(TBH_table) <- fullnames


# Make tables at each phylogenetic level
tribe_table <- combine_otus("Tribe", TBH_table, TBH_taxonomy)
clade_table <- combine_otus("Clade", TBH_table, TBH_taxonomy)
lineage_table <- combine_otus("Lineage", TBH_table, TBH_taxonomy)
order_table <- combine_otus("Order", TBH_table, TBH_taxonomy)
class_table <- combine_otus("Class", TBH_table, TBH_taxonomy)
phylum_table <- combine_otus("Phylum", TBH_table, TBH_taxonomy)


# Combine tables at each level into one giant table
full_table <- rbind(TBH_table, tribe_table, clade_table, lineage_table, order_table, class_table, phylum_table)

# Remove groups that are unclassified at any level
classified <- grep("unclassified", rownames(full_table))
classified1 <- grep("__$", rownames(full_table))
input_table <- full_table[-c(classified, classified1),]
input_table <- t(input_table)
input_table <- as.data.frame(input_table)

#Group samples based on time from spring mixing. Date of spring mixing in Trout Bog obtained from the NSIDC Global Lake and River Phenology Database

TBH_dates <- extract_date(rownames(input_table))
TBH_years <- substr(rownames(input_table), start = 9, stop = 10)

daynum <- c()
for(i in 1:length(TBH_dates)){
  if (TBH_years[i]  == "05"){
    daynum[i] <- as.numeric(TBH_dates[i] - extract_date(c("TBH14Apr05")))
  }else if (TBH_years[i] == "07"){
    daynum[i] <- as.numeric(TBH_dates[i] - extract_date(c("TBH13Apr07")))
  }else if (TBH_years[i] == "08"){
    daynum[i] <- as.numeric(TBH_dates[i] - extract_date(c("TBH30Apr08")))
  }else if (TBH_years[i] == "09"){
    daynum[i] <- as.numeric(TBH_dates[i] - extract_date(c("TBH21Apr09")))
  }
}

#Remove winter samples

input_table <- input_table[which(daynum > 0),]
daynum <- daynum[which(daynum > 0)]

#Place each daynumber into categories based on date from spring mixing - every fortnight?

seasons <- seq(from = 0, to = 220, by = 40)

period <- c()
for(i in 1:length(daynum)){
  period[i] <- length(which(seasons <= daynum[i]))
}
period <- factor(period, levels = 1:6)

#Now that I have my groups and input table, I'm ready to run the indicator analysis

time.indicators <- multipatt(x = input_table, cluster = period, func = "r.g", control = how(nperm = 9999))

# Save data
write.csv(time.indicators$sign, file = "C:/Users/amlin/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/indicators_by_times.csv")
