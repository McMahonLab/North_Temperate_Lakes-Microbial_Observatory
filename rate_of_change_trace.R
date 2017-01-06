#Rate of change over imte

library(OTUtable)
library(reshape2)
#library(dplyr)
library(ggplot2)
library(grid)
library(raster)
#library(ggrepel)
library(ape)
library(phyloseq)
library(vegan)

seq_table <- read.csv("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/bog_seqstable.csv", row.names = 1)
seq_taxonomy <- read.csv(file = "C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/seqs.98.cleantaxonomy.csv", row.names = 1, header = T, colClasses = c("character"))

seqs <- read.dna("C:/Users/Alex/Desktop/North_Temperate_Lakes-Microbial_Observatory/Data/16S_data/qc.bogs.clean.min25.fasta", format = "fasta")
d <- dist.dna(seqs, model = "raw")
bogtree <- nj(d)

# Set up phyloseq object to run UniFrac on
seq_table <- remove_reps(seq_table)
OTU <- otu_table(as.matrix(seq_table), taxa_are_rows = T)
TAX <- tax_table(as.matrix(seq_taxonomy))
sampledata <- sample_data(data.frame(Bog = substr(colnames(seq_table), start = 1, stop = 2), Layer = substr(colnames(seq_table), start = 3, stop = 3), Year = substr(colnames(seq_table), start = 9, stop = 10), row.names = colnames(seq_table), stringsAsfactors = F)) 
alldata <- phyloseq(OTU, TAX, sampledata, bogtree)

lakes <- c("FB", "CB", "WS", "NS", "TB", "SS", "MA", "HK")
layers <- c("E", "H")
years <- c("05", "07", "08", "09")

for(a in 1:length(lakes)){
  for(b in 1:length(layers)){
    for(c in 1:length(years)){
      # Only perform the analysis of the lake-layer-year subset actually exists
      if(length(which(sampledata$Bog == lakes[a] & sampledata$Layer == layers[b] & sampledata$Year == years[c])) > 0){
        subset <- prune_samples(sampledata$Bog == lakes[a] & sampledata$Layer == layers[b] & sampledata$Year == years[c], alldata)
        # Calculate UniFrac and convert to long form
        x <- UniFrac(subset, weighted = T, normalize = T)
        x1 <- melt(as.matrix(x))
        # Make a list of unique dates so that I can select only distances between consecutive samples from the long UniFrac table
        samples <- unique(x1$Var1)
        samples <- samples[order(extract_date(samples))]
        x2 <- x1[1, ]
        for(i in 1:(length(samples)-1)){
          new.x <- x1[which(x1$Var1 == samples[i] & x1$Var2 == samples[i+1]), ]
          x2 <- rbind(x2, new.x)
        }
        x2 <- x2[2:length(samples), ]
        # Use the ice-off date as time zero
        if(lakes[a] == "CB"){
          if(years[c] == "05"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("CBE12APR05")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("CBE12APR05")))
          }else if(years[c] == "07"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("CBE12APR07")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("CBE12APR07")))
          }else if(years[c] == "08"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("CBE25APR08")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("CBE25APR08")))
          }else if(years[c] == "09"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("CBE20APR09")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("CBE20APR09")))
          }
        }else if(lakes[a] == "TB"){
          if(years[c] == "05"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("TBE15APR05")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("TBE15APR05")))
          }else if(years[c] == "07"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("TBE18APR07")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("TBE18APR07")))
          }else if(years[c] == "08"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("TBE30APR08")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("TBE30APR08")))
          }else if(years[c] == "09"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("TBE22APR09")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("TBE22APR09")))
          }
        }else{
          if(years[c] == "05"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("CBE13APR05")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("CBE13APR05")))
          }else if(years[c] == "07"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("CBE16APR07")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("CBE16APR07")))
          }else if(years[c] == "08"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("CBE27APR08")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("CBE27APR08")))
          }else if(years[c] == "09"){
            x2$Adj.Date1 <- as.numeric(extract_date(x2$Var1) - extract_date(c("CBE21APR09")))
            x2$Adj.Date2 <- as.numeric(extract_date(x2$Var2) - extract_date(c("CBE21APR09")))
          }
        }
        # Add a column of time difference to my new table
        x2$Time <- x2$Adj.Date2 - x2$Adj.Date1
        
        # Rate of change analysis
        x2$Rate <- x2$value/x2$Time
        x2 <- x2[which(x2$Rate != Inf),]
        assign(paste("rate_",lakes[a], layers[b], years[c] , sep = ""), x2)
        
        
      }
    }
  }
}

rate_CBE07$Year <- factor(substr(rate_CBE07$Var1, start = 9, stop = 10), levels = c("07", "09"))
rate_CBE09$Year <- factor(substr(rate_CBE09$Var1, start = 9, stop = 10), levels = c("07", "09"))
CBE_rates <- rbind(rate_CBE07[, c(5, 7, 8)], rate_CBE09[, c(5, 7, 8)])
ggplot(data = CBE_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "CBE")

rate_CBH07$Year <- factor(substr(rate_CBH07$Var1, start = 9, stop = 10), levels = c("07", "09"))
rate_CBH09$Year <- factor(substr(rate_CBH09$Var1, start = 9, stop = 10), levels = c("07", "09"))
CBH_rates <- rbind(rate_CBH07[, c(5, 7, 8)], rate_CBH09[, c(5, 7, 8)])
ggplot(data = CBH_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "CBH")

rate_NSE07$Year <- factor(substr(rate_NSE07$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
rate_NSE08$Year <- factor(substr(rate_NSE08$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
rate_NSE09$Year <- factor(substr(rate_NSE09$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
NSE_rates <- rbind(rate_NSE07[, c(5, 7, 8)], rate_NSE08[, c(5, 7, 8)], rate_NSE09[, c(5, 7, 8)])
ggplot(data = NSE_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "NSE")

rate_NSH07$Year <- factor(substr(rate_NSH07$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
rate_NSH08$Year <- factor(substr(rate_NSH08$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
rate_NSH09$Year <- factor(substr(rate_NSH09$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
NSH_rates <- rbind(rate_NSH07[, c(5, 7, 8)], rate_NSH08[, c(5, 7, 8)], rate_NSH09[, c(5, 7, 8)])
ggplot(data = NSH_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "NSH")

rate_SSE07$Year <- factor(substr(rate_SSE07$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
rate_SSE08$Year <- factor(substr(rate_SSE08$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
rate_SSE09$Year <- factor(substr(rate_SSE09$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
SSE_rates <- rbind(rate_SSE07[, c(5, 7, 8)], rate_SSE08[, c(5, 7, 8)], rate_SSE09[, c(5, 7, 8)])
ggplot(data = SSE_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "SSE")

rate_SSH07$Year <- factor(substr(rate_SSH07$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
rate_SSH08$Year <- factor(substr(rate_SSH08$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
rate_SSH09$Year <- factor(substr(rate_SSH09$Var1, start = 9, stop = 10), levels = c("07", "08", "09"))
SSH_rates <- rbind(rate_SSH07[, c(5, 7, 8)], rate_SSH08[, c(5, 7, 8)], rate_SSH09[, c(5, 7, 8)])
ggplot(data = SSH_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "SSH")

rate_TBE05$Year <- factor(substr(rate_TBE05$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_TBE07$Year <- factor(substr(rate_TBE07$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_TBE08$Year <- factor(substr(rate_TBE08$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_TBE09$Year <- factor(substr(rate_TBE09$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
TBE_rates <- rbind(rate_TBE05[, c(5, 7, 8)], rate_TBE07[, c(5, 7, 8)], rate_TBE08[, c(5, 7, 8)], rate_TBE09[, c(5, 7, 8)])
ggplot(data = TBE_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "TBE")

rate_TBH05$Year <- factor(substr(rate_TBH05$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_TBH07$Year <- factor(substr(rate_TBH07$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_TBH08$Year <- factor(substr(rate_TBH08$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_TBH09$Year <- factor(substr(rate_TBH09$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
TBH_rates <- rbind(rate_TBH05[, c(5, 7, 8)], rate_TBH07[, c(5, 7, 8)], rate_TBH08[, c(5, 7, 8)], rate_TBH09[, c(5, 7, 8)])
ggplot(data = TBH_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "TBH")

rate_MAE05$Year <- factor(substr(rate_MAE05$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_MAE07$Year <- factor(substr(rate_MAE07$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_MAE08$Year <- factor(substr(rate_MAE08$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_MAE09$Year <- factor(substr(rate_MAE09$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
MAE_rates <- rbind(rate_MAE05[, c(5, 7, 8)], rate_MAE07[, c(5, 7, 8)], rate_MAE08[, c(5, 7, 8)], rate_MAE09[, c(5, 7, 8)])
ggplot(data = MAE_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "MAE")

rate_MAH05$Year <- factor(substr(rate_MAH05$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_MAH07$Year <- factor(substr(rate_MAH07$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_MAH08$Year <- factor(substr(rate_MAH08$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
rate_MAH09$Year <- factor(substr(rate_MAH09$Var1, start = 9, stop = 10), levels = c("05", "07", "08", "09"))
MAH_rates <- rbind(rate_MAH05[, c(5, 7, 8)], rate_MAH07[, c(5, 7, 8)], rate_MAH08[, c(5, 7, 8)], rate_MAH09[, c(5, 7, 8)])
ggplot(data = MAH_rates, aes(x = Adj.Date2, y = Rate, group = Year, color = Year)) + geom_path() + theme_bw() + labs(title = "MAH")

