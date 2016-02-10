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
clean_taxonomy <- function(taxonomy_file, table){
  y <- read.table(taxonomy_file, header = T)
  keep <- match(rownames(table), y$OTU, nomatch = NA)
  y <- y[keep,]
  return(y)
}

# Convert taxonomy file from a single string per OTU to having each taxonomic level in a separate column
expand_taxa <- function(taxonomy_file, remove_bootstrap){
  tax <- read.csv(file = taxonomy_file, header = T, row.names = 1, colClasses = c("character"))
  
  sp <- strsplit(tax$Taxonomy, split = ";")
  kingdom <- c()
  phyla <- c()
  class <- c()
  order <- c()
  lineage <- c()
  clade <- c()
  tribe <- c()
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
  if(remove_bootstrap == T){
    kingdom <- gsub("\\(\\d*\\)", "", kingdom)
    phyla <- gsub("\\(\\d*\\)", "", phyla)
    class <- gsub("\\(\\d*\\)", "", class)
    order <- gsub("\\(\\d*\\)", "", order)
    lineage <- gsub("\\(\\d*\\)", "", lineage)
    clade <- gsub("\\(\\d*\\)", "", clade)
    tribe <- gsub("\\(\\d*\\)", "", tribe)
  }
  
  new.tax <- data.frame(kingdom, phyla, class, order, lineage, clade, tribe, stringsAsFactors=F)
  rownames(new.tax) <- rownames(tax)
  colnames(new.tax) <- c("Kingdom", "Phylum", "Class", "Order", "Lineage", "Clade", "Tribe")
  return(new.tax)
}

