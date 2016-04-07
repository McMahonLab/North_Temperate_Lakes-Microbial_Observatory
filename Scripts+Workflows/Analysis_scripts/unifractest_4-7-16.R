library(OTUtable)
library(phyloseq)
library(ggplot2)
library(ape)
library(vegan)
library(grid)
library(reshape2)
path2repo <- "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Figures/"
data(otu_table)
data(taxonomy)
data(metadata)

seqs <- read.dna("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_repseqs_07Jul15.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

# Make OTU table, taxonomy, and sampledata datasets
OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(taxonomy))

sampledata <- sample_data(data.frame(Bog = substr(colnames(otu_table), start = 1, stop = 2), Layer = substr(colnames(otu_table), start = 3, stop = 3), Year = substr(colnames(otu_table), start = 9, stop = 10), row.names = colnames(otu_table), stringsAsfactors = F))                                                                                                                                          

alldata <- phyloseq(OTU, TAX, sampledata, bogtree)
years <- c("05", "07", "08", "09")

hypodata <- prune_samples(sampledata$Layer == "H", alldata)

x <- UniFrac(hypodata, weighted = T, normalize = T)
x <- as.matrix(x)

MAsamples <- x[which(substr(rownames(x), start = 1, stop = 2) != "MA"),which(substr(rownames(x), start = 1, stop = 2) == "MA")]
mean.Unifrac <- 1 - rowSums(MAsamples)/dim(MAsamples)[2]
sampledates <- extract_date(names(mean.Unifrac))
samplelakes <- substr(names(mean.Unifrac), start = 1, stop = 2)
uf <- data.frame(mean.Unifrac, sampledates, samplelakes)
colnames(uf) <- c("meanUF", "Date", "Lake")

#Visually test each lake
ggplot(data = uf[which(uf$Lake == "TB"),], aes(x = Date, y = meanUF)) + geom_line() #Good
ggplot(data = uf[which(uf$Lake == "CB"),], aes(x = Date, y = meanUF)) + geom_line() #No trend
ggplot(data = uf[which(uf$Lake == "NS"),], aes(x = Date, y = meanUF)) + geom_line() #Unclear, mixing event visible
ggplot(data = uf[which(uf$Lake == "SS"),], aes(x = Date, y = meanUF)) + geom_line() #Probably in some years not all
ggplot(data = uf[which(uf$Lake == "HK"),], aes(x = Date, y = meanUF)) + geom_line() #Outliers drive (-) trend
ggplot(data = uf[which(uf$Lake == "FB"),], aes(x = Date, y = meanUF)) + geom_line() #No trend, some outliers
ggplot(data = uf[which(uf$Lake == "WS"),], aes(x = Date, y = meanUF)) + geom_line() #No trend

#Add Bray-Curtis to the mix

bc <- distance(hypodata, "bray")
bc <- as.matrix(bc)
MAsamples <- bc[which(substr(rownames(bc), start = 1, stop = 2) != "MA"),which(substr(rownames(bc), start = 1, stop = 2) == "MA")]
meanBC <- 1 - rowSums(MAsamples)/dim(MAsamples)[2]
uf$meanBC <- meanBC

ggplot(data = uf[which(uf$Lake == "TB"),], aes(x = Date, y = meanBC)) + geom_line() #really unclear
ggplot(data = uf[which(uf$Lake == "CB"),], aes(x = Date, y = meanBC)) + geom_line() #trend in 07, not 09
ggplot(data = uf[which(uf$Lake == "NS"),], aes(x = Date, y = meanBC)) + geom_line() # (-) trends
ggplot(data = uf[which(uf$Lake == "SS"),], aes(x = Date, y = meanBC)) + geom_line() #some trend w outliers
ggplot(data = uf[which(uf$Lake == "HK"),], aes(x = Date, y = meanBC)) + geom_line() #  (-) trend
ggplot(data = uf[which(uf$Lake == "FB"),], aes(x = Date, y = meanBC)) + geom_line() #(-) trend
ggplot(data = uf[which(uf$Lake == "WS"),], aes(x = Date, y = meanBC)) + geom_line() #No trend

ggplot(data = uf, aes(x = meanBC, y = meanUF)) + geom_point()
#Not great agreement. What is different now? Lack of phylogenetic info in Bray-Curtis at OTU level?

#Try calculating clade based Bray-Curtis

clade_table <- combine_otus("Clade", otu_table, taxonomy)
hypo_clade_table <- bog_subset("..H", clade_table)
bc.clade <- vegdist(t(hypo_clade_table), method = "bray")
bc.clade <- as.matrix(bc.clade)
bc.clade <- bc.clade[which(substr(rownames(bc.clade), start = 1, stop = 2) != "MA"), which(substr(rownames(bc.clade), start = 1, stop = 2) == "MA" )]

uf$cladeBC <- 1 - rowSums(bc.clade)/dim(bc.clade)[2]

ggplot(data = uf, aes(x = meanBC, y = cladeBC)) + geom_point() #Good
ggplot(data = uf, aes(x = meanUF, y = cladeBC)) + geom_point() #Not as good

#Do the evenness outliers play into this?

even <- estimate_richness(hypodata, measures = "Shannon")/log(length(which(taxa_sums(hypodata) > 0)))
even <- even[which(substr(rownames(even), start = 1, stop = 2) != "MA"),]
uf$even <- even

ggplot(data = uf, aes(x = meanUF, y = cladeBC, color = even)) + geom_point() #Low evenness at low end for both metrics
ggplot(data = uf, aes(x = meanUF, y = cladeBC, color = Lake)) + geom_point() 
ggplot(data = uf, aes(x = meanUF, y = cladeBC, color = Date)) + geom_point()

#Nothing obvious here. Try adding in the representative Mary samples and seeing how those agree with the clade Bc

MAHtable <- bog_subset("MAH", otu_table)
MAHtable05 <- year_subset("05", MAHtable)
MAHtable07 <- year_subset("07", MAHtable)
MAHtable08 <- year_subset("08", MAHtable)
MAHtable09 <- year_subset("09", MAHtable)

otu_table2 <- otu_table
otu_table2$MAH01REP05 <- rowSums(MAHtable05)/dim(MAHtable05)[2]
otu_table2$MAH01REP07 <- rowSums(MAHtable07)/dim(MAHtable07)[2]
otu_table2$MAH01REP08 <- rowSums(MAHtable08)/dim(MAHtable08)[2]
otu_table2$MAH01REP09 <- rowSums(MAHtable09)/dim(MAHtable09)[2]


OTU2 <- otu_table(as.matrix(otu_table2), taxa_are_rows = T)
sampledata2 <- sample_data(data.frame(Bog = substr(colnames(otu_table2), start = 1, stop = 2), Layer = substr(colnames(otu_table2), start = 3, stop = 3), Year = substr(colnames(otu_table2), start = 9, stop = 10), row.names = colnames(otu_table2), stringsAsfactors = F)) 
alldata_reps <- phyloseq(OTU2, TAX, sampledata2, bogtree)

hypodata_reps <- prune_samples(sampledata2$Layer == "H", alldata_reps)
x2 <- UniFrac(hypodata_reps, weighted = T, normalize = T)
x2 <- as.matrix(x2)

sim <- c()
for(i in 1:690){
  sampleyear <- substr(rownames(x2)[i], start = 9, stop = 10)
  if (sampleyear == "05"){
    sim[i] <- x2[i, 691]
  }else if (sampleyear == "07"){
    sim[i] <- x2[i, 692]
  }else if (sampleyear == "08"){
    sim[i] <- x2[i, 693]
  }else if (sampleyear == "09"){
    sim[i] <- x2[i, 694]
  }
}
sim <- 1 - sim[which(substr(rownames(x2), start = 1, stop = 2) != "MA")]

uf$repUF <- sim

ggplot(data = uf, aes(x = meanUF, y = repUF)) + geom_point() #Very nice
ggplot(data = uf, aes(x = cladeBC, y = repUF)) + geom_point() #Not as nice, but not bad either

ggplot(data = uf, aes(x = even, y = repUF, color = Lake)) + geom_point()
ggplot(data = uf, aes(x = even, y = cladeBC, color = Lake)) + geom_point()
ggplot(data = uf, aes(x = Date, y = repUF, color = even)) + geom_point(size = 3)

uf$Year <- substr(rownames(uf), start = 9, stop = 10)
ggplot(data = uf[which(uf$Lake == "TB" & uf$Year == "07"),], aes(x = Date, y = repUF, color = even)) + geom_point(size = 3)

# A pretty strict evenness threshold instead of a mixing threshold may fix a lot of my problems. 
ggplot(data = uf[which(uf$Lake == "TB" & uf$Year == "07"),], aes(x = Date, y = even)) + geom_point(size = 3)
ggplot(data = uf[which(uf$Lake == "TB" & uf$Year == "07"),], aes(x = even, y = repUF)) + geom_point(size = 3)
# Mixed samples have even < 0.4 in TB. What about NS?

ggplot(data = uf[which(uf$Lake == "NS" & uf$Year == "08"),], aes(x = Date)) + geom_line(aes(y = repUF)) + geom_line(aes(y = cladeBC), color = "blue")
# Too strict here - more like 0.35

ggplot(data = uf[which(uf$Lake == "HK" & uf$Year == "07"),], aes(x = Date, y = even)) + geom_point(size = 3)
ggplot(data = uf[which(uf$Lake == "HK" & uf$Year == "07"),], aes(x = even, y = repUF)) + geom_point(size = 3)
# 0.55 would work here

ggplot(data = uf[which(uf$Lake == "NS" & uf$Year == "07"),], aes(x = cladeBC, y = meanUF, color = even)) + geom_point(size = 3)
ggplot(data = uf[which(uf$Lake == "CB" & uf$Year == "07"),], aes(x = even, y = repUF)) + geom_point(size = 3)

#Switch my rep UF to weighted = F and try again
#Clade BC and mean BC match up perfectly with my old data (visually). The UF metrics do not. 

x3 <- UniFrac(hypodata_reps, weighted = F, normalize = T)
x3 <- as.matrix(x3)

sim <- c()
for(i in 1:690){
  sampleyear <- substr(rownames(x3)[i], start = 9, stop = 10)
  if (sampleyear == "05"){
    sim[i] <- x3[i, 691]
  }else if (sampleyear == "07"){
    sim[i] <- x3[i, 692]
  }else if (sampleyear == "08"){
    sim[i] <- x3[i, 693]
  }else if (sampleyear == "09"){
    sim[i] <- x3[i, 694]
  }
}
sim <- 1 - sim[which(substr(rownames(x3), start = 1, stop = 2) != "MA")]

uf$repUF_unweighted <- sim

ggplot(data=uf, aes(y = repUF, x = repUF_unweighted)) + geom_point()
ggplot(data=uf, aes(y = repUF, x = cladeBC)) + geom_point()
ggplot(data=uf, aes(y = repUF_unweighted, x = cladeBC)) + geom_point()

#Unweighted looks much more linear in comparison to cladeBC

ggplot(data = uf[which(uf$Lake == "HK" & uf$Year == "07"),], aes(x = Date)) + geom_line(aes(y = repUF_unweighted)) + geom_line(aes(y = cladeBC), color = "blue") 

#Are the outliers in HK also uneven?
ggplot(data = uf[which(uf$Lake == "HK" & uf$Year == "07"),], aes(x = Date)) + geom_point(aes(y = repUF_unweighted, color = even))
hist(uf$even[which(uf$Lake == "HK" & uf$Year == "07")])
hist(uf$even[which(uf$Lake == "NS" & uf$Year == "07")])
hist(uf$even[which(uf$Lake == "SS" & uf$Year == "07")])
hist(uf$even[which(uf$Lake == "TB" & uf$Year == "07")])

#Yes

#New plan: calculated significance in the linear trends using unweighed representative sample UniFrac and mean cladeBC, while excluding mixed dates and points with evenness less than 0.3

#Sidenote: I'm  not sure I previously calculated evenness quite right. 

even <- estimate_richness(hypodata, measures = "Shannon")/log(estimate_richness(hypodata, measures = "Observed"))
even <- even[which(substr(rownames(even), start = 1, stop = 2) != "MA"),]
uf$even <- even

#Good correlation, looks like I was just a little off. Set threshold at 0.5 now

# add mixing variables to my dataframe

cont.mixes <- c()
for(i in 1:dim(uf)[1]){
  find <- grep(substr(rownames(uf)[i], start = 1, stop = 10), metadata$Sample_Name)
  if(length(find) > 0){
    samples <- metadata[find,]
    top <- samples$DO[which(samples$Depth == "0.5")]
    if(length(top) == 0){
      top <- samples$DO[which(samples$Depth == "1")]
    }
    bottom <- min(samples$DO)
    cont.mixes[i] <- top - bottom
  }
}

uf$Mixing.cont <- cont.mixes
uf$Mixing.bin <- uf$Mixing.cont < 5.5

uf$Year <- substr(rownames(uf), start = 9, stop = 10)
#Julian Date variables

julian <- c()
for(i in 1:dim(uf)[1]){
  if(uf$Year[i] == "05"){
    julian[i] <- as.numeric(uf$Date[i] - extract_date(c("TBH01Jun05")))
  }else if(uf$Year[i] == "07"){
    julian[i] <- as.numeric(uf$Date[i] - extract_date(c("TBH01Jun07")))
  }else if(uf$Year[i] == "08"){
    julian[i] <- as.numeric(uf$Date[i] - extract_date(c("TBH01Jun08")))
  }else if(uf$Year[i] == "09"){
    julian[i] <- as.numeric(uf$Date[i] - extract_date(c("TBH01Jun09")))
  }
}

uf$DayNum <- julian

#Go for the home run first. Samples must have DO diff > 4, evenness > 0.5, and DayNum less than 153 (to exclude November onwards)
#allow variation by lake and year
#I expect that dimictic lakes will have a significant positive correlation between unweighted representative sample UniFrac and DayNum, and other lakes will not (unless they acted like a dimictic lake that year)

model1 <- lm(repUF_unweighted ~  -1 + Lake:DayNum:Year:Mixing.bin, data=uf[which(uf$even > 0.5 & uf$DayNum < 153),])
summary(model1)
#I'm still seeing a variety of trends in the 2007 lakes especially. Plot these out.

ggplot(data = uf[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$Lake == "CB"& uf$Year == "07"& uf$DayNum < 153),], aes(x = DayNum, y = repUF_unweighted)) + geom_line() #(+) but model says (-) - are the september time points throwing it off? DO diff of 5.2
ggplot(data = uf[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$Lake == "FB"& uf$Year == "07"& uf$DayNum < 153),], aes(x = DayNum, y = repUF_unweighted)) + geom_line() #(+) but model says (-)
ggplot(data = uf[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$Lake == "HK"& uf$Year == "07"& uf$DayNum < 153),], aes(x = DayNum, y = repUF_unweighted)) + geom_line() #(-) but model says (+) - still one outlier
ggplot(data = uf[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$Lake == "NS"& uf$Year == "07"& uf$DayNum < 153),], aes(x = DayNum, y = repUF_unweighted)) + geom_line() #no trend, model syas (-)
ggplot(data = uf[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$Lake == "SS"& uf$Year == "07"& uf$DayNum < 153),], aes(x = DayNum, y = repUF_unweighted)) + geom_line() #no trend, model says (+)
ggplot(data = uf[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$Lake == "TB" & uf$Year == "07"& uf$DayNum < 153),], aes(x = DayNum, y = repUF_unweighted)) + geom_line() #(+), model says (+) - mixed samples still included?
ggplot(data = uf[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$Lake == "WS"& uf$Year == "07"& uf$DayNum < 153),], aes(x = DayNum, y = repUF_unweighted)) + geom_line() #slight (+)? model says no trend


# Compare to clade BC
model2 <- lm(cladeBC ~  Lake:DayNum:Year , data=uf[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153),])
summary(model2)

#Do I get different results if I manually compare lakes?

#HK 07
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "HK"& uf$Year == "07")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "HK"& uf$Year == "07")]
summary(lm( v1 ~ v2)) # p = 0.0593
        
#FB 07
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "FB"& uf$Year == "07")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "FB"& uf$Year == "07")]
summary(lm( v1 ~ v2)) # p = 0.999

#CB 07
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "CB"& uf$Year == "07")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "CB"& uf$Year == "07")]
summary(lm( v1 ~ v2)) # p = 0.000114

#WS 07
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "WS"& uf$Year == "07")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "WS"& uf$Year == "07")]
summary(lm( v1 ~ v2)) # p = 0.038189

#NS 07
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "NS"& uf$Year == "07")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "NS"& uf$Year == "07")]
summary(lm( v1 ~ v2)) # p = 3e-5

#TB 07
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "TB"& uf$Year == "07")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "TB"& uf$Year == "07")]
summary(lm( v1 ~ v2)) # p = 7e-6

#SS 07
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "SS"& uf$Year == "07")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "SS"& uf$Year == "07")]
summary(lm( v1 ~ v2)) # p = 0.22

#SS 08
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "SS"& uf$Year == "08")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "SS"& uf$Year == "08")]
summary(lm( v1 ~ v2)) # p = 0.22

#SS 09
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "SS"& uf$Year == "09")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "SS"& uf$Year == "09")]
summary(lm( v1 ~ v2)) # p = 0.22

#TB 05
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "TB"& uf$Year == "05")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "TB"& uf$Year == "05")]
summary(lm( v1 ~ v2)) # p = 0.083

#TB 08
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "TB"& uf$Year == "08")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "TB"& uf$Year == "08")]
summary(lm( v1 ~ v2)) # p = 3e-5

#TB 09
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "TB"& uf$Year == "09")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "TB"& uf$Year == "09")]
summary(lm( v1 ~ v2)) # p = 0.587

#CB 09
v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "CB"& uf$Year == "09")]
v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$even > 0.5 & uf$DayNum < 153 & uf$Lake == "CB"& uf$Year == "09")]
summary(lm( v1 ~ v2)) # p = 0.000114

#Write a function to take bog and year and measure its trend. I want to know sample size, r^2, p-value, and slope

phenology <- function(bog, year){
  v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$Lake == bog & uf$Year == year)]
  v2 <- uf$DayNum[which(uf$Mixing.bin == F & uf$Lake == bog & uf$Year == year)]
  model <- lm(v1 ~ v2)
  r <- cor(v2, v1)
  n <- length(v2)
  p <- summary(model)$coefficients[2, 4]
  m <- model$coefficients[2]
  b <- model$coefficients[1]
  range <- range(v1)[2] - range(v1)[1]
print("Correlation: "); print(r)
print("Sample Size: "); print(n)
print("p-value: "); print(p)
print("slope: "); print(m)
print("intercept: "); print(b)
print("range:"); print(range)
plot(v2, v1)
}

v1 <- uf$repUF_unweighted[which(uf$Mixing.bin == F & uf$Lake == bog & uf$Year == year)]
v2 <- uf$meanBC[which(uf$Mixing.bin == F & uf$Lake == bog & uf$Year == year)]
plot(v2, v1)
