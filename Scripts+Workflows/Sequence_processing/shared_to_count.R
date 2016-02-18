cat("Input shared file:")
shared <- scan(what=character(), nmax=1)
cat("Output count file:")
out.count <- scan(what=character(), nmax=1)

x <- read.delim(file=shared, header=T)
table <- x[,4:length(x[1,])]
header <- substr(colnames(table), start=2, stop=7)
real <- which(header > 0)
header <- header[real]

table <- t(table[,real])
total <- rowSums(table)
columns <- c("Representative_Sequence", "total", as.character(x$Group))
count_table <- data.frame(header, total, table, row.names=NULL)
colnames(count_table) <- columns
write.table(count_table, file=out.count, quote=F, sep="\t", row.names=F, col.names=T)
