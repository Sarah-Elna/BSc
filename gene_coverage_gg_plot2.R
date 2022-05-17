# Gene cover gg plot in R, but written by Sarah E. Valentin

require(gplots)
install.packages("heatmap.plus")

setwd("C:/Users/Sarah/Documents/AU/6/Bachelor/GenomeDK/")

sample.filename = "seq_lengths.txt"

sample.data= as.matrix(read.table(sample.filename,header=T,row.names=1))
sample.len = sample.data[2:nrow(sample.data),]
reference.len = as.numeric(sample.data[1,])