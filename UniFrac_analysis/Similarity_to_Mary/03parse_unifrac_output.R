#Parse unifrac output into something more useful
args <- commandArgs(trailingOnly=TRUE)
#Subset samples wanted
library(OTUtable)
library(reshape2)

unifrac <- read.delim(paste("C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/UniFrac_analysis/", args[1], "-", args[2], "-", args[3], ".trewsummary", sep = ""), sep = "\t", header=T, check.names=F)

#split up groups column into sample 1 and sample 2

split.groups <- strsplit(as.character(unifrac$Groups), "-")
split.groups <- unlist(split.groups)
sample1 <- split.groups[seq(1, by=2, len=dim(unifrac)[1])]
sample2 <- split.groups[seq(2, by=2, len=dim(unifrac)[1])]

unifrac2 <- data.frame(sample1, sample2, unifrac$WScore)
#I only want sample comparisons between different bogs.
bog1 <- substr(unifrac2$sample1, start=1, stop=2)
bog2 <- substr(unifrac2$sample2, start=1, stop=2)
unifrac2 <- unifrac2[which(bog1 != bog2),]
#Now I want to average by trout bog sample
library(dplyr)
x <- group_by(unifrac2, sample2)
product <- summarize(x, mean(unifrac.WScore))
product$date <- extract_date(product$sample2)
colnames(product) <- c("sample", "UnifracW", "Date")
product$UnifracW <- 1-product$UnifracW

#Output results - I'll use these for plotting
write.csv(product, paste(args[1], "-", args[2], "-", args[3], "-UniFrac.csv", sep=""))

#Take samples from stratified dates only
data(metadata)
metalakes <- substr(metadata$Sample_Name, start = 1, stop = 3)
metayears <- substr(metadata$Sample_Name, start = 9, stop = 10)

metabog <- metabog <- metadata[which(metalakes == args[1] & metayears == args[3]), c(1,2,4)]
metabog2 <- dcast(metabog, Sample_Name~Depth, fun.aggregate=mean)
temp.range <- c()
for(i in 1:dim(metabog2)[1]){
  sample <- as.numeric(metabog2[i,])
  temp.range[i] <- sample[2] - min(sample[which(is.na(sample) == F)])
}
mix.dates <- extract_date(metabog2$Sample_Name[which(temp.range < 1)])
remove <- match(mix.dates, product$Date)

stratified <- product[!product$Date %in% mix.dates,]
#quantify phenology
n <- dim(stratified)[1]
time <- as.numeric(stratified$Date)
r <- cor(stratified$UnifracW, time)
model <- lm(stratified$UnifracW~time)
p <- summary(model)$coefficients[2, 4]

print("Samples: "); print(n)
print("r^2:  "); print(r)
print("p-value: "); print(p)
