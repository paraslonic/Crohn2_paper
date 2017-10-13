library("reshape2")
library("data.table")

t = fread("recombinations_recent.txt",skip=1)
strains = fread("strains", col.names = c("num","strain"))

t. = subset(t, t$RecipientStrain == 2)
plot(t.$End, t.$DonorLineage+1, type = "n", ylim=c(0,5))
for(i in 1:nrow(t.)){
  lines(c(t.[i,]$Start, t.[i,]$End), 
        c(t.[i,]$DonorLineage,t.[i,]$DonorLineage), 
        col = t.[i,]$DonorLineage, lwd = 4)
  
}

t. = subset(t, t$RecipientStrain == 1)
for(i in 1:nrow(t.)){
  lines(c(t.[i,]$Start, t.[i,]$End), 
        c(t.[i,]$DonorLineage+0.2,t.[i,]$DonorLineage+0.2), 
        col = t.[i,]$DonorLineage, lwd = 4)
}


t.[, sum(End - Start), by=DonorLineage]

## global
t
t[,length:=End-Start]
t[,sum(length),by=.(StrainName,DonorLineage)]
hist(t$`log(BF)`, xlim = c(0,500), breaks  =1000)
t.ok = t[`log(BF)` > 200]
t.ok[,sum(length),by=.(RecipientStrain,DonorLineage)]
View(t)
#
plot(t$length, t$`log(BF)`, log="xy")
rec = t[, sum(End - Start), by=.(RecipientStrain,DonorLineage)]
recm = dcast(rec, RecipientStrain ~ DonorLineage)
recm[is.na(recm)] = 0
library("gplots")
heatmap.2(as.matrix(recm[,-1]), scale = "row")
