library("IRanges")
library("GenomicRanges")
library("reshape2")
library("gplots")

setwd("blastout")
blastFiles <- list.files(".",".blast")
names = gsub(".blast","",blastFiles)

f <- lapply(blastFiles, read.table)
col_names = c("qid", "sid", "identity", "length", "mismatches", "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "querycoveragepersubject", "subjectlength")
f <- lapply(f, function(x) {colnames(x) <- col_names; x})

ldata <-  list()
for (i in 1:length(f)){
  df <- as.data.frame(f[[i]])
  plasmid_length <- data.frame(plasmid=df$sid, length=df$subjectlength)
  plasmid_length <- plasmid_length[!duplicated(plasmid_length$plasmid),]
  for (j in 1:nrow(df)){
    if (df$sstart[j] > df$send[j]){
      sstart <- df$sstart[j]
      df$sstart[j] <- df$send[j]
      df$send[j] <- sstart
    }
  }
  plranges <- GRanges(seqnames = df$sid, ranges = IRanges(start = df$sstart, end = df$send))
  plranges <- reduce(plranges)
  t <- data.frame(id = seqnames(plranges), width = width(plranges))
  t <- aggregate(t$width, by=list(t$id), sum)
  colnames(t) <- c("plasmid","coverage")
  t <- merge(t, plasmid_length, by="plasmid", all.x = TRUE)
  t$coverage <- t$coverage/t$length
  t$sample <- names[i]
  ldata[[i]] <- t
}

table <- do.call(rbind, ldata)
wtable <- dcast(table, plasmid ~ sample, value.var = "coverage")
wtable[is.na(wtable)] <- 0
rownames(wtable) <- wtable$plasmid
pdf("../plasmid_coverage_heatmap.pdf")
heatmap.2(t(as.matrix(wtable[,-1])), margins = c(6,6), trace="none", col=colorRampPalette(c("gray30", "dodgerblue2")))
dev.off()
write.table(wtable, "../plasmid_coverage_table.txt",quote = F,sep="\t")
