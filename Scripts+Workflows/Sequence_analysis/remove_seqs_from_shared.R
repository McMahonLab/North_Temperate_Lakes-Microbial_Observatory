cat("Input shared file:")
shared <- scan(what=character(), nmax=1)
cat("Input names file:")
names <- scan(what=character(), nmax=1)

cat("Output shared file:")
out.shared <- scan(what=character(), nmax=1)

x <- read.delim(file=shared, header=T)
table <- x[,4:length(x[1,])]
header <- substr(colnames(table), start=2, stop=7)
colnames(table) <- header

y <- read.delim(file=names, header=F)
realseqs <- y$V1

keep <- match(realseqs, header, nomatch=NA)
label <- x[,1]
numOtus <- rep(length(keep), length(label))
Group <- x[,2]
new_names <- c("label", "Group", "numOtus", header[keep])
new_shared <- data.frame(label, Group, numOtus, table[,keep])
colnames(new_shared) <- new_names

write.table(new_shared, file=out.shared, quote=F, sep="\t", row.names=F, col.names=T)

