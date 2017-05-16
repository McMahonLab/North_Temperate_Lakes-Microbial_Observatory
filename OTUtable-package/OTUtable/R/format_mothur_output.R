# Clean up shared file into a traditional OTU table

# 1. remove label, numOTus
# 2. transpose
# 3. shorten sample names

clean_shared <- function(shared_file, trim.names){
  x <- read.table(shared_file, header = T, row.names = 2)
  x$label <- NULL
  x$numOtus <- NULL
  if(trim.names == T){
    split.names <- strsplit(rownames(x), split = "\\.")
    samples <- c()
    for(i in 1:length(split.names)){
      samples[i] <- split.names[[i]][1]
    }
    rownames(x) <- make.unique(samples)
  }
  x <- t(x)
  return(x)
}

# Keep only the OTU column of the taxonomy file, and only OTUs that survived the subsampling step in the OTU table
clean_mothur_taxonomy <- function(taxonomy_file, table, remove_bootstrap){
  
  tax <- read.table(file=taxonomy_file, header=T, row.names=1, colClasses=c("character"))
  tax$Size <- NULL
  sp <- strsplit(tax$Taxonomy, split = ";")
  kingdom <- c(rep("", length(tax$Taxonomy)))
  phyla <- c(rep("", length(tax$Taxonomy)))
  class <- c(rep("", length(tax$Taxonomy)))
  order <- c(rep("", length(tax$Taxonomy)))
  lineage <- c(rep("", length(tax$Taxonomy)))
  clade <- c(rep("", length(tax$Taxonomy)))
  tribe <- c(rep("", length(tax$Taxonomy)))
                    
  for(i in 1:length(sp)){
     group <- sp[i][[1]]
     kingdom[i] <- group[1]
     phyla[i] <- group[2]
     class[i] <- group[3]
     order[i] <- group[4]
     lineage[i] <- group[5]
     clade[i] <- group[6] 
     tribe[i] <- group[7]
     }
                    
  y <- data.frame(kingdom, phyla, class, order, lineage, clade, tribe, stringsAsFactors=F) 
  rownames(y) <- rownames(tax)
  colnames(y) <- c("Kingdom", "Phylum", "Class", "Order", "Lineage", "Clade", "Tribe") 
  
  keep <- match(rownames(y), rownames(table))
  y <- y[which(is.na(keep) == F), ]
  y <- y[order(rownames(y)), ]
  
  if(remove_bootstrap == T){
    y$Kingdom <- gsub("\\(\\d*\\)", "", y$Kingdom)
    y$Phylum <- gsub("\\(\\d*\\)", "", y$Phylum)
    y$Class <- gsub("\\(\\d*\\)", "", y$Class)
    y$Order <- gsub("\\(\\d*\\)", "", y$Order)
    y$Lineage <- gsub("\\(\\d*\\)", "", y$Lineage)
    y$Clade <- gsub("\\(\\d*\\)", "", y$Clade)
    y$Tribe <- gsub("\\(\\d*\\)", "", y$Tribe)
  } 
  return(y)
}

clean_TaxAss_taxonomy <- function(taxonomy_file, table, remove_bootstrap){
  tax <- read.csv(taxonomy_file, header = T, row.names = 1)
  colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Lineage", "Clade", "Tribe")
  keep <- match(rownames(tax), rownames(table))
  tax <- tax[which(is.na(keep) == F), ]
  y <- tax[order(rownames(tax)), ]
  
  if(remove_bootstrap == T){
    y$Kingdom <- gsub("\\(\\d*\\)", "", y$Kingdom)
    y$Phylum <- gsub("\\(\\d*\\)", "", y$Phylum)
    y$Class <- gsub("\\(\\d*\\)", "", y$Class)
    y$Order <- gsub("\\(\\d*\\)", "", y$Order)
    y$Lineage <- gsub("\\(\\d*\\)", "", y$Lineage)
    y$Clade <- gsub("\\(\\d*\\)", "", y$Clade)
    y$Tribe <- gsub("\\(\\d*\\)", "", y$Tribe)
  } 
  return(y)
}
                  