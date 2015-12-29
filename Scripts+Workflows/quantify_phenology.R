#Quantify the observed phenology in community similarity between dimictic and meromictic hypolimnia

#Set up the environment
library(reshape2)
library(OTUtable)
library(vegan)
data(otu_table)
data(taxonomy)
data(metadata)
clade_table <- combine_otus("Clade", otu_table, taxonomy)
mary <- bog_subset("MAH", clade_table)

#Specify lakes and years for testing
dimictic <- c("NSH", "TBH", "SSH")
years <- c("05", "07", "08", "09")

#I will remove mixed dates from the analysis. In order to do this, I first need to identify mixed dates
#Use these vectors to search metadata for mixed dates in the for loop
metalakes <- substr(metadata$Sample_Name, start=1, stop=3)
metayears <- substr(metadata$Sample_Name, start=9, stop=10)

#Set up empty matrices to hold data
#m = slope
#b = intercept
#r = r^2 correlation coefficient
#p = p-value
#n = sample size
m <- matrix(NA, length(dimictic), length(years))
r <- matrix(NA, length(dimictic), length(years))
n <- matrix(NA, length(dimictic), length(years))
b <- matrix(NA, length(dimictic), length(years))
p <- matrix(NA, length(dimictic), length(years))

#Repeat this analysis for every dimictic lake/year
for(i in 1:length(dimictic)){
  #Get dimictic lake
  bog <- bog_subset(dimictic[i], clade_table)
  
  for(j in 1:length(years)){
    #Get year of data
    mary.year <- year_subset(years[j], mary)
    data <- year_subset(years[j], bog)
    #Get year's metadata
    metabog <- metadata[which(metalakes == dimictic[i] & metayears == years[j]), c(1,2,4)]
    #If statement necessary because there is no data for NSH and SSH in 2005
    if(dim(data)[2] > 0){
      #Find dates with less than 1 degree of difference between 0.5m and the max sampling depth (defined as mixed)
      metabog2 <- dcast(metabog, Sample_Name~Depth, fun.aggregate=mean)
      mixes <- extract_date(metabog2$Sample_Name[which(metabog2$"0.5" - metabog2[,dim(metabog2)[2]-1] < 1)])
      bogdates <- extract_date(colnames(data))
      #Remove these dates from dataset
      parsed.data <- data[,which(is.na(match(bogdates, mixes)) == T)]
      parsed.dates <- extract_date(colnames(parsed.data))
      parsed.data <- parsed.data[,order(parsed.dates)]
      final.dates <- extract_date(colnames(parsed.data))
      #Calculate Bray-Curtis Similarity and average by dimictic sample
      mean.bc <- c()
      for(k in 1:dim(parsed.data)[2]){
        bc.dis <- c()
        for(l in 1:dim(mary.year)[2]){
          test <- cbind(parsed.data[,k], mary.year[,l])
          bc.dis[l] <- 1 - vegdist(t(test), method="bray")
        }
        mean.bc[k] <- mean(bc.dis) 
      }
      #Fill relevant information about the trend into the appropriate place in each matrix
      time <- as.numeric(final.dates - final.dates[1])
      r[i,j] <- cor(time, mean.bc)
      n[i,j] <- length(mean.bc)  
      bestfit <- lm(mean.bc ~ time)
      m[i,j] <- bestfit$coefficients[2]
      b[i,j] <- bestfit$coefficients[1]
      p[i,j] <- summary(bestfit)$coefficients[2,4]
    }
  }
}
#Label rows and columnx in the matrices
colnames(m) <- years
colnames(r) <- years
colnames(p) <- years
colnames(b) <- years
rownames(m) <- dimictic
rownames(r) <- dimictic
rownames(p) <- dimictic
rownames(b) <- dimictic
