dimictic <- c("NSH", "TBH", "SSH")
years <- c("05", "07", "08", "09")

library(reshape2)
library(OTUtable)
library(vegan)
data(otu_table)
data(taxonomy)
data(metadata)
clade_table <- combine_otus("Clade", otu_table, taxonomy)
metalakes <- substr(metadata$Sample_Name, start = 1, stop = 3)
metayears <- substr(metadata$Sample_Name, start = 9, stop = 10)
mary <- bog_subset("MAH", clade_table)

m <- matrix(NA, length(dimictic), length(years))
r <- matrix(NA, length(dimictic), length(years))
n <- matrix(NA, length(dimictic), length(years))
b <- matrix(NA, length(dimictic), length(years))
p <- matrix(NA, length(dimictic), length(years))
for(i in 1:length(dimictic)){
  bog <- bog_subset(dimictic[i], clade_table)
  
  for(j in 1:length(years)){
    mary.year <- year_subset(years[j], mary)
    data <- year_subset(years[j], bog)
    metabog <- metadata[which(metalakes == dimictic[i] & metayears == years[j]), c(1,2,4)]
    if(dim(data)[2] > 0){
      metabog2 <- dcast(metabog, Sample_Name~Depth, fun.aggregate=mean)
      mixes <- extract_date(metabog2$Sample_Name[which(metabog2$"0.5" - metabog2[,dim(metabog2)[2]-1] < 1)])
      bogdates <- extract_date(colnames(data))
      parsed.data <- data[,which(is.na(match(bogdates, mixes)) == T)]
      parsed.dates <- extract_date(colnames(parsed.data))
      parsed.data <- parsed.data[,order(parsed.dates)]
      final.dates <- extract_date(colnames(parsed.data))
      mean.bc <- c()
      for(k in 1:dim(parsed.data)[2]){
        bc.dis <- c()
        for(l in 1:dim(mary.year)[2]){
          test <- cbind(parsed.data[,k], mary.year[, l])
          bc.dis[l] <- 1 - vegdist(t(test), method="bray")
        }
        mean.bc[k] <- mean(bc.dis) 
      }
      time <- as.numeric(final.dates - final.dates[1])
      r[i, j] <- cor(time, mean.bc)
      n[i, j] <- length(mean.bc)  
      bestfit <- lm(mean.bc ~ time)
      m[i, j] <- bestfit$coefficients[2]
      b[i, j] <- bestfit$coefficients[1]
      p[i, j] <- summary(bestfit)$coefficients[2, 4]
    }
  }
}
colnames(m) <- years
colnames(r) <- years
colnames(p) <- years
colnames(b) <- years
rownames(m) <- dimictic
rownames(r) <- dimictic
rownames(p) <- dimictic
rownames(b) <- dimictic
