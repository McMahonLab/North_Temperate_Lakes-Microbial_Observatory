cat("Input file name:")
filein <- scan(what=character(), nmax=1)
cat("Output file name:")
fileout <- scan(what=character(), nmax=1)
x <- read.delim(file=filein, header=T, row.names=2)
table <- x[,3:length(x[1,])]
header <- substr(colnames(table), start=2, stop=7)

names <- data.frame(header, header)
write.table(names, file=fileout, quote=F, sep="\t", row.names=F, col.names=F)
