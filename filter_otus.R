library(OTUtable)
data(otu_table)

filter_otus <- function(table, abundance, persistence){
  # Set up objects for persistance data
  pers <- c()
  num_samples <- dim(table)[2]
  
  # Set up objects for abundance data
  meets_abun_threshold <- c()
  abun_threshold <- sum(table[, 1]) * (abundance/100)
  
  # Loop through each OTU
  for(i in 1:dim(table)[1]){
    # Calculate the percentage of samples present for each OTU
    present <- length(which(table[i, ] > 0))
    pers[i] <- present/num_samples * 100
    
    # Calculate whether or not each OTU crosses the relative abundance threshold
    meets_abun_threshold[i] <- length(which(table[i, ] >= abun_threshold)) > 0
  }
  
  # Keep only OTUs that meet the set requirements
  filtered_table <- table[which(pers >= persistence & meets_abun_threshold == T), ]
  return(filtered_table)
}

#Example
library(OTUtable)
data(otu_table)
new_table <- filter_otus(otu_table, abundance = 0.1, persistence = 50)
