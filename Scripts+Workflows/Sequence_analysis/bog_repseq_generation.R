#Goal: make an OTU fasta file

#I need: a key that says what seqIDs went into which OTU
#something that says what the representative sequence for each OTU is
#the seqID fasta file

#After some digging, the taxonomy reported was a consensus of all seqs in an OTU. No rep seqs were chosen

#I only need the list file (which seqs in which otu) and the seqID fasta

list <- read.table("D:/server files/EMB_bog_tags/deblurring_w_r/qc.bogs.clean.min25.phylip.an.list", header=T, row.names=1, fill=T, colClasses = c("character"))

fasta <- read.table("D:/server files/EMB_bog_tags/deblurring_w_r/qc.bogs.clean.min25.fasta", header=F, colClasses = c("character"))


#Approach:
#- loop through 0.02 row
#- split up mulitple entries and choose one randomly if necessary
#- once seq id is chosen, add > and match to entry in fasta file
#- get sequence and replace seqID with otu number
#- create and export a fasta file with otus and rep seqs

#Remove OTUs that didn't make it through subsampling - about 100
#Load taxonomy file to see which OTUs to keep

taxonomy <- read.csv("C:/Users/Alex/Dropbox/Deblurred bog tags/bogs_OTU_taxonomy07Jan15.csv", header=T)
otus <- as.character(taxonomy[,1])
keep <- match(otus, colnames(list))
cand <- list[3, keep]

new_fasta <-c()
for(i in 1:length(cand)){
  entry <- cand[i]
  if(nchar(entry) > 6){
    all <- strsplit(as.character(entry), split = ",")
    entry <- sample(all[[1]], 1)
  }
  search <- paste(">", entry, sep = "")
  seqid <- match(search, fasta$V1)
 
  new_fasta <- append(new_fasta, paste(">", colnames(list)[keep[i]], sep=""), length(new_fasta))
  new_fasta <- append(new_fasta, fasta$V1[seqid+1], length(new_fasta))
}

write.table(new_fasta, file="bog_repseqs_07Jul15.fasta", quote=F, row.names=F, col.names=F)



