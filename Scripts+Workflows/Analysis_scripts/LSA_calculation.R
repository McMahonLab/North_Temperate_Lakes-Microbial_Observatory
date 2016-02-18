setwd("C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/Additional_analyses/Networks")
minreads <- 100
minlsa <- .75

library(OTUtable)
data(otu_table)
data(taxonomy)

x <- otu_table
#  Optional year subset
#  x <- year_subset("08", x)
#  Optional month subset
#  x <- x[,which(substr(colnames(x), start=6, stop=8) == "MAY")]
#  Optional bog subset
x <- bog_subset("MAH", x)
x <- x[which(rowSums(x) >= minreads), ]

lsa_prep <- function(table){
  table <- table[,which(substr(colnames(table), start=12, stop=13) != "R2")]
  table <- table[which(rowSums(table)>0), ]
  dates <- as.Date(substr(colnames(table), start=4, stop=10), format="%d%b%y")
  table <- table[,order(dates)]
  ztable <- matrix(NA, ncol=length(table[1, ]), nrow=length(table[, 1]))
  for (i in 1:length(table[, 1])){
    otu <- as.numeric(table[i, ])
    zotu <- (otu- mean(otu))/sd(otu)
    ztable[i,] <- zotu
  }
  return(ztable)
}

temp <- lsa_prep(x)
write.table(temp, file="temp.txt", row.names=F, col.names=F, sep = "\t")
out <- "MAH_network.txt"

# My desktop path
# system(paste("D:/fastlsa_win/fastlsa.exe -i temp.txt", " -o ", out, " -d 0", " -m ", minlsa, sep=""))
# My laptop path
system(paste("C:/Users/amlinz16/fastlsa_win/fastlsa.exe -i temp.txt", " -o ", out, " -d 0", " -m ", minlsa, sep=""))

network <- read.table(file=out, header=T)
network$index1 <- rownames(x)[network$index1 + 1]
network$index2 <- rownames(x)[network$index2 + 1]
write.table(network, file=out, row.names=F, col.names=T, sep = "\t", quote=F)

# Number edges:
dim(network)[1]
# Number nodes:
length(unique(c(network$index1, network$index2)))
