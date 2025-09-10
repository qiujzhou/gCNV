library(rCNV)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(parallel)
library(ggpubr)

##### parallelization and basic files
num_cores = 20
cl <- makeCluster(num_cores)
probe2gene <- readRDS("path2file") #### two columns table: 1st col: probe name,2nd:gene name
gff <- read.table("path2file",h=T) #### info table of probes: chr, start, end, probe name
info <- read.table("path2info",h=T) #### info table of samples: ID,pop,lat,long,cluster,etc.


############ pca on non-CNV snps
cline.evec <- read.table("path2PCA_result.evec") #### PCA using plink
cline.evec <- cbind(cline.evec,info[match(cline.evec$id,info$indi),])
cline.evec <- cline.evec[!is.na(cline.evec$cluster),]
cline.evec <- cline.evec[match(rownames(depth_mean60.norm.scale),cline.evec$id),]
cline.eval <- read.table("path2PCA_result.eval")
cline.eval$V1[1:2]/sum(cline.eval$V1)

# color by genetic cluster
p5.cluster <- ggplot(cline.evec) + geom_point(aes(PC1,PC2,col=cluster))+xlab("PC1 (6.1%)") +
  scale_color_manual(values= c(brewer.pal(9,"Set1"))[c(2:6,9,7,8,1)]) + 
  theme_classic()+
  labs(color="Genetic cluster") + ylab("PC2 (2.2%)") +
  theme(
    panel.grid.major = element_blank(),    # Remove major gridlines
    panel.grid.minor = element_blank() 
  ) +
  ggtitle("Single-copy SNPs") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))
p5.cluster

# color by longitude
p5 <- ggplot(cline.evec) + geom_point(aes(-V2,-V3,col=Long))+xlab("PC1 (6.1%)") +
  labs(color="Longitude") + ylab("PC2 (2.2%)") +
  theme_classic()+
  ggtitle("Single-copy SNPs") +
  scale_color_continuous(low = "purple", high = "orange") +
  theme(
    panel.grid.major = element_blank(),    # Remove major gridlines
    panel.grid.minor = element_blank() 
  ) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))
p5



############################### normalization ##################
depth_mean60 <- readRDS("path2depth") #### mean depth for each probe each sample,
                                      #### averaged across 60bp in the center of each probe
depth_mean60[is.na(depth_mean60)] <- 0
depth_mean60 <- depth_mean60[,intersect(info$indi,colnames(depth_mean60))]
hist(colSums(depth_mean60))
hist(rowSums(depth_mean60),breaks = 100)
rm(depth_mean60)

### remove probe too low average depth
sum(rowSums(depth_mean60)<ncol(depth_mean60)*5)
depth_mean60 <- depth_mean60[rowSums(depth_mean60) > ncol(depth_mean60)*5,]
dim(depth_mean60)

### MedR normalization account for library size variation
pseudo <- apply(depth_mean60, 1, function(xx) {
  exp(mean(log(as.numeric(xx)[as.numeric(xx) > 0])))
})
nf <- apply(depth_mean60, 2, function(xx) {
  median(as.numeric(xx)/pseudo, na.rm = T)
})

depth_mean60.norm <- round(depth_mean60/rep(nf, each = dim(depth_mean60)[1]), 2)
depth_mean60.norm[is.na(depth_mean60.norm)] <- 0

### PCA on all probes, filter out outliers 
depth_mean60.norm <- t(depth_mean60.norm)
depth_mean60.norm.scale <- scale(depth_mean60.norm)
depth_mean60.pca <- prcomp(depth_mean60.norm.scale,scale. = T,center = T,retx = T)
depth_mean60.rotation <- cbind(depth_mean60.pca$x,info[match(rownames(depth_mean60.norm),info$indi),])

# plot PCA
ggplot(depth_mean60.rotation) + geom_point(aes(PC1,PC2,col=cluster))+
  scale_color_manual(values= c(brewer.pal(9,"Set1"))[c(2:6,9,7,8,1)])
# remove outlier
depth_mean60.norm <- depth_mean60.norm[rownames(depth_mean60.rotation)[depth_mean60.rotation$PC1 > -100],] 
depth_mean60.norm.scale <- depth_mean60.norm.scale[rownames(depth_mean60.rotation)[depth_mean60.rotation$PC1 > -100],] # remove outlier
depth_mean60 <- depth_mean60[,rownames(depth_mean60.rotation)[depth_mean60.rotation$PC1 > -100]] # remove outlier
dim(depth_mean60.norm.scale)

### add stats to info table
info <- info[match(rownames(depth_mean60.norm),info$indi),]
info$lib <- rowSums(depth_mean60.norm)
info$var <- apply(depth_mean60.norm,1,function(x){sd(x)/mean(x)})
info$proj <- unlist(lapply(info$indi,function(x){unlist(strsplit(x,"_P"))[1]})) ### sequencing batch


############################################################ PCA based normalization
depth_mean60.svd <- La.svd(as.matrix(depth_mean60.norm.scale))
depth_mean60.svd$d[1:10]^2/sum(depth_mean60.svd$d^2)

### check correlation between PCs and confounding factors
pc.con.cor <- lapply(1:5,function(x){
  return(data.frame(PC=x,
                    proj = summary(lm(depth_mean60.svd$u[,x] ~ info$proj))$adj.r.squared,
                    lib.size = summary(lm(depth_mean60.svd$u[,x] ~ info$lib))$adj.r.squared,
                    doc.var = summary(lm(depth_mean60.svd$u[,x] ~ info$var))$adj.r.squared,
                    pop.stru = summary(lm(depth_mean60.svd$u[,x] ~ cline.evec$PC1 + cline.evec$PC2))$adj.r.squared,
                    probe.effi = summary(lm(depth_mean60.svd$vt[x,] ~ colMeans(depth_mean60.norm)))$adj.r.squared,
                    probe.gc = summary(lm(depth_mean60.svd$vt[x,] ~ 
                                            probe.stat$X5_pct_gc[match(colnames(depth_mean60.norm),probe.stat$probe)]))$adj.r.squared))
})
pc.con.cor <- do.call(rbind,pc.con.cor)

### PCA based normalization, remove PC1 and PC5
depth_mean60.svd$d[c(1,5)] <- 0 
depth_mean60.svd.norm <- depth_mean60.svd$u %*% diag(depth_mean60.svd$d) %*% depth_mean60.svd$vt
dimnames(depth_mean60.svd.norm) <- dimnames(depth_mean60.norm)
saveRDS(depth_mean60.svd.norm,"depth_mean60.norm.rds")


################################################################### extract CNV probes
cnv.dep60.norm <- depth_mean60.svd.norm[,intersect(cnv_probe,colnames(depth_mean60.svd.norm))]
probe_mean <- colMeans(depth_mean60.norm[rownames(cnv.dep60.norm),])
probe_sd <- apply(depth_mean60.norm[rownames(cnv.dep60.norm),],2,function(x){sd(x)})
depth_final <- apply(depth_mean60.svd.norm, 1, function(x){
  x*probe_sd + probe_mean
})  #### unscale normalized depth table

saveRDS(depth_final,"depth_final.rds")

