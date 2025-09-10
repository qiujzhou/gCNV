library(rCNV)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(parallel)
library(ggpubr)
library(mclust)

##### parallelization and basic files
num_cores = 20
cl <- makeCluster(num_cores)
probe2gene <- readRDS("path2file") #### two columns table: 1st col: probe name,2nd:gene name
gff <- read.table("path2file",h=T) #### info table of probes: chr, start, end, probe name
info <- read.table("path2info",h=T) #### info table of samples: ID,pop,lat,long,cluster,etc.



##################################################################### Mclust filtering
cnv.mclust <- parApply(cnv.dep.final,1,function(x){
  library(mclust)
  hc <- hc(x, modelName = "E", use = "STD")
  BIC1 <- mclustBIC(x, initialization = list(hcPairs = hc))
  hc <- hc(x, modelName = "V", use = "STD")
  BIC3 <- mclustBIC(x, initialization = list(hcPairs = hc))
  BIC <- mclustBICupdate(BIC1, BIC3)
  
  best <- Mclust(x,x = BIC)
  y <- tapply(x, best$classification, mean)
  #names(y) <- sort(unique(best$classification))
  return(y[as.factor(best$classification)])
},cl = cl)

rownames(cnv.mclust) <- colnames(cnv.dep.final)

#max.copy <- apply(cnv.mclust, 2, function(x){length(unique(x))})
#barplot(table(max.copy),main = "P.abies-P.obovata")

##### filter based on mclust result
mclust.filter <- apply(cnv.mclust, 2, function(x){
  y <- table(x)
  if (length(y) == 1 ){
    FALSE
  } else if(length(y) == 2 & min(y) < 5){
    FALSE
  } else {
    TRUE
  }
})
sum(mclust.filter)

cnv.dep.final <- cnv.dep.final[mclust.filter,]
saveRDS(cnv.dep.final,"/home/zhouqj/proj5000/cline/cnv.dep.final.rds")
cnv.probe.final <- rownames(cnv.dep.final)

###### nonCNV genes
nonCNV.gene <- unique(probe2gene$gids[!probe2gene$gids %in% unique(probe2gene$gids[probe2gene$probe %in% cnv.probe.final])])
nonCNV.probe <- probe2gene$probe[probe2gene$gids %in% nonCNV.gene]
nonCNV.gene <- unique(probe2gene$gids[probe2gene$probe %in% nonCNV.probe])
###### CNV genes
cnv.gene <- unique(probe2gene$gids[probe2gene$probe %in% cnv.probe.final])

################################################################# SNP filtering
snp.pos <- fread("filepath") #### chr,pos
snp.pos$probe <- parApply(snp.pos,1,function(x, probe ){
  if (any(probe$start < as.numeric( x[2]) & probe$end > as.numeric(x[2]) & probe$chr == as.character(x[1])) ) {
    return(probe$probe[probe$start < as.numeric( x[2]) & probe$end > as.numeric(x[2]) & probe$chr == as.character(x[1])][1])
  } else{
    return(FALSE)
  }
},cl = cl,probe=gff)

snp.pos$gene <- probe2gene$gids[match(snp.pos$probe,probe2gene$probe)]
snp.singcopy.pos <- snp.pos[snp.pos$gene %in% nonCNV.gene,]  #### SNPs in single-copy region
snp.cnv.pos <- snp.pos[snp.pos$gene %in% cnv.gene,]   #### SNPs in gCNV
write.table(snp.singcopy.pos[,1:2],"snp.singcopy.pos", 
            quote = F, col.names = F,row.names = F)
write.table(snp.cnv.pos[,1:2],"snp.cnv.pos", 
            quote = F, col.names = F,row.names = F)

####### final gCNV table
cnv.dep.final <- depth_final[intersect(cnv.probe.final,rownames(depth_final)),]
saveRDS(cnv.dep.final,"cnv.dep.final.rds")


