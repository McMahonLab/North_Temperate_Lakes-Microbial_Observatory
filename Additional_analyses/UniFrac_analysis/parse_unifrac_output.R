# Parse unifrac output

unifrac <- read.delim("C:/Users/amlinz16/Desktop/North_Temperate_Lakes-Microbial_Observatory/Additional_analyses/UniFrac_analysis/temp.trewsummary", sep = "\t", header=T, check.names=F)

# split up groups column into sample 1 and sample 2

split.groups <- strsplit(as.character(unifrac$Groups), "-")
split.groups <- unlist(split.groups)
sample1 <- split.groups[seq(1, by=2, len=dim(unifrac)[1])]
sample2 <- split.groups[seq(2, by=2, len=dim(unifrac)[1])]

unifrac2 <- data.frame(sample1, sample2, unifrac$WScore)
# I only want sample comparisons between different bogs.
bog1 <- substr(unifrac2$sample1, start=1, stop=2)
bog2 <- substr(unifrac2$sample2, start=1, stop=2)
unifrac2 <- unifrac2[which(bog1 != bog2),]
# Now I want to average by trout bog sample
library(dplyr)
x <- group_by(unifrac2, sample2)
product <- summarize(x, mean(unifrac.WScore))
product$date <- extract_date(product$sample2)
colnames(product) <- c("sample", "UnifracW", "Date")
product$UnifracW <- 1-product$UnifracW
library(ggplot2)
ggplot(data=product, aes(x=Date, y=UnifracW)) + geom_line() + coord_cartesian(xlim=extract_date(c("TBE01Jun07", "TBE01Dec07")))

# Same trend observed as in Sorenson Index - increase over time, then decrease during mixing event