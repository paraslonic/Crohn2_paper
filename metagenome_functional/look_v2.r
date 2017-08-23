library("modeest")


setwd("/mnt/yakov/crohn2/text/metagenome_functional")

pa1 <- read.delim("kazan5500_2014_07_10_4_L01_GorK_F3_subsetAsmet2_pathabundance.tsv", head = T)
colnames(pa1) <- c("pathway","ab1")

pa2 <- read.delim("kazan5500_2016_12_14_1_L06_FHM_F3_pathabundance.tsv", head = T)
colnames(pa2) <- c("pathway","ab2")

pc1 <- read.delim("kazan5500_2014_07_10_4_L01_GorK_F3_subsetAsmet2_pathcoverage.tsv")
colnames(pc1) <- c("pathway","cov1")

pc2 <- read.delim("kazan5500_2016_12_14_1_L06_FHM_F3_pathcoverage.tsv")
colnames(pc2) <- c("pathway","cov2")

t <- merge(pa1, pa2,   all = TRUE)
t <- subset(t, !grepl("UN",t$pathway))
t <- subset(t, !grepl("__",t$pathway))

# determine pseudo count value as maximum nonzero value
q <- t[,-1]
q[is.na(q)] <- max(q,na.rm = T)
pseudo_count <- min(q)
t[is.na(t)] <- pseudo_count

# scale normalize 

ratio <- t$ab1/t$ab2
mode <- mlv(ratio)
t$ab1 <- t$ab1/mode[1]$M

t <- merge(t, pc1, all.x=TRUE)
t <- merge(t, pc2, all.x=TRUE)

t <- t[order(t$ab1/t$ab2, decreasing = T),]
t[is.na(t)] <- 0
t <- subset(t, apply(t[,c(4,5)],1,max)>0.1)

write.table(t, "table.txt", quote = F, sep="\t")

View(t)
