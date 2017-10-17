library("reshape2")
library("data.table")
library("gplots")

t = fread("recombinations_recent.txt",skip=1)
strains = fread("strains", col.names = c("num","strain"))
lineages = fread("lineages", col.names = c("id","lineage"))

t[,length:=End-Start]
t[,sum(length),by=.(StrainName,DonorLineage)]
rec = as.data.frame(dcast(t, RecipientStrain ~ DonorLineage, fun.aggregate = sum, value.var = "length"))
rec = rec[,-1]

colnames(rec) = lineages[as.integer(colnames(rec))]$lineage
rownames(rec) = strains[as.integer(rownames(rec))]$strain
heatmap.2(as.matrix(rec), trace="none", 
          col = colorRampPalette(c("black","dark orange","orange","yellow","white")),
          margins=c(14,14), cexCol = 1.2, cexRow=1, Colv=FALSE, Rowv=FALSE, dendrogram = "none")


ani = read.delim("oani.csv")
rownames(ani) = ani$name
ani = ani[,-1]
as.numeric(ani)
ani[is.na(ani)] = 0
ani = ani+t(ani)
diag(ani) = 100
ani = apply(ani, c(1,2),as.numeric)
ani = 100 - ani
hc = as.dendrogram(hclust(as.dist(ani)))

clade_order <- order.dendrogram(hc)
clade_name <- labels(hc)
clade_position <- data.frame(clade_name,clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name,row.names(rec))
rec = rec[new_order,]

