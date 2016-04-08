library(OTUtable)
library(vegan)
library(phyloseq)

data(otu_table)
data(taxonomy)
data(metadata)

seqs <- read.dna("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_repseqs_07Jul15.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

MAHtable <- bog_subset("MAH", otu_table)
MAHtable05 <- year_subset("05", MAHtable)
MAHtable07 <- year_subset("07", MAHtable)
MAHtable08 <- year_subset("08", MAHtable)
MAHtable09 <- year_subset("09", MAHtable)


REP05 <- rowSums(MAHtable05)/dim(MAHtable05)[2]
REP07 <- rowSums(MAHtable07)/dim(MAHtable07)[2]
REP08 <- rowSums(MAHtable08)/dim(MAHtable08)[2]
REP09 <- rowSums(MAHtable09)/dim(MAHtable09)[2]

#Just test the 2007 sample to cut down on run time

test <- seq(from = 250, to = 999, by = 15)
UFW.typeI <- c()
UF.typeI <- c()
BC.typeI <- c()
UFW.typeII <- c()
UF.typeII <- c()
BC.typeII <- c()

B = 100

for(k in 1:B){
  
  fakesamples07 <- matrix(ncol = 50, nrow = 6208)
  for(i in 1:length(test)){
    fakesample <- c(rep(0, 6208))
    matched.otus <- sample(which(REP07 > 0), size = test[i],  replace = F)
    fakesample[matched.otus] <- REP07[matched.otus]
    remaining.counts <- 5000 - sum(fakesample)
    other.otus <- sample(rownames(otu_table), size = remaining.counts, prob = rowSums(otu_table), replace = T)
    fake.counts <- table(other.otus)
    fakesample[match(names(fake.counts), names(REP07))] <- fakesample[match(names(fake.counts), names(REP07))] + fake.counts
    fakesample[which(is.na(fake.sample) == T)] <- 0
    fakesamples07[,i] <- fakesample
  }
  fakesamples07 <- data.frame(fakesamples07, REP07)
  rownames(fakesamples07) <- rownames(otu_table)
  
  # Make phyloseq object
  OTU <- otu_table(as.matrix(fakesamples07), taxa_are_rows = T)
  TAX <- tax_table(as.matrix(taxonomy))
  testdata07 <- phyloseq(OTU, TAX, bogtree)
  
  
  # Type 1 error - how often does each metric produce p > 0.5 when there is a trend?
  
  x <- UniFrac(testdata07, weighted = T, normalize = T)
  x <- as.matrix(x)
  UFW07 <- x[1:50,51]
  UFW.typeI[k] <- summary(lm(UFW07 ~ test))$coefficients[2, 4]
  
  y <- UniFrac(testdata07, weighted = F, normalize = T)
  y <- as.matrix(y)
  UF07 <- y[1:50,51]
  UF.typeI[k] <- summary(lm(UF07 ~ test))$coefficients[2, 4]
  
  z <- vegdist(t(fakesamples07), method = "bray")
  z <- as.matrix(z)
  BC07 <- z[1:50,51]
  BC.typeI[k] <- summary(lm(BC07 ~ test))$coefficients[2, 4]
  
  #How often does type II error - detection of a trend when there is none - occur?

  model <- lm(UFW07 ~ sample(1:50, size = 50, replace = F))
  UFW.typeII[k] <- summary(model)$coefficients[2,4]
  
  model <- lm(UF07 ~ sample(1:50, size = 50, replace = F))
  UF.typeII[k] <- summary(model)$coefficients[2,4]
  
  model <- lm(BC07 ~ sample(1:50, size = 50, replace = F))
  BC.typeII[k] <- summary(model)$coefficients[2,4]
  print(k)
}

print("UniFrac, weighted: type I")
print(length(which(UFW.typeI >= 0.05))/B)

print("UniFrac, unweighted: type I")
print(length(which(UF.typeI >= 0.05))/B)

print("Bray Curtis: type I")
print(length(which(BC.typeI >= 0.05))/B)

print("UnIIFrac, weighted: type II")
print(length(which(UFW.typeII < 0.05))/B)

print("UnIIFrac, unweighted: type II")
print(length(which(UF.typeII < 0.05))/B)

print("Bray CurtIIs: type II")
print(length(which(BC.typeII < 0.05))/B)

UFW.typeII <- c()
UF.typeII <- c()
BC.typeII <- c()

for(k in 1:1000){
  shuffle <- sample(1:50, size = 50, replace = F)
  model <- lm(UFW07 ~ shuffle)
  UFW.typeII[k] <- summary(model)$coefficients[2,4]
  
  model <- lm(UF07 ~ shuffle)
  UF.typeII[k] <- summary(model)$coefficients[2,4]
  
  model <- lm(BC07 ~ shuffle)
  BC.typeII[k] <- summary(model)$coefficients[2,4]
}


