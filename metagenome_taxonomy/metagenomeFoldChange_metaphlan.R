library(ggplot2)
library(ggrepel)

setwd("/mnt/yakov/crohn2/metagenome/metaphlan/")

p1 = read.delim("s1.metaphlan2_out", head = F,skip = 1)
p2 = read.delim("s2.metaphlan2_out", head = F,skip = 1)

p1 = subset(p1, grepl("s__",p1$V1) & !grepl("t__",p1$V1))
p2 = subset(p2, grepl("s__",p2$V1) & !grepl("t__",p2$V1))

p1$V1 = gsub(".*s__","",p1$V1)
p2$V1 = gsub(".*s__","",p2$V1)

t = merge(p1,p2, by = "V1", all = T)
colnames(t) = c("name","abundance1","abundance2")

t[is.na(t)] = min(c(p1$V2, p2$V2))/10 # pseudo count
#t$abundance1 = t$abundance1/sum(t$abundance1)
#t$abundance2 = t$abundance2/sum(t$abundance2)

# select with 
tmax = apply(t[,c(2,3)], 1, max)
t = t[tmax > 0.1,]

# plot with fold change
fold_change = t$abundance1/t$abundance2
mean_abundance = apply(t[,c(2,3)], 1, mean)
t$fold_change = fold_change
t$mean_abundance = mean_abundance
t = subset(t, !grepl("Homo",t$name))
t = subset(t, !grepl("phage",t$name))
#t = subset(t, !grepl("phage",t$name))
t$name= gsub('\\[',"",t$name)
t$name = gsub("\\]","",t$name)

makeShortName = function(x){
  p = strsplit(x," ")[[1]]
  if(length(p) != 2 || grepl("uncultured",p[1])) return(x)
  paste0(substr(p[1],1,1), ". ", p[2])
}
t$name = sapply(t$name, makeShortName)

ggplot(data = t, aes(x = fold_change, y = mean_abundance)) + theme_classic() + 
  geom_text_repel(aes(label = name),size = 3, box.padding = unit(0.3, "lines"), segment.alpha = 0.3) +
  geom_point(colour = "dodgerblue4", size = 3)+scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='log2') + labs(x = "abundance fold change", y="mean abundance")
  
### table
t = t[order(t$fold_change, decreasing = TRUE),]
write.table(t, "metagenome_abundance_fold_change_metaphlan.txt",quote = F, sep = "\t",row.names=F)

