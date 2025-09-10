library(rCNV)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(parallel)
library(ggpubr)
library(geodist)
library(raster)
library(rgdal)
library(terra)


######################################### Fst
fst = read.table("all_pop_fst.txt",header = F)      
colnames(fst) <- c("pop1","pop2","fst")
poplist = unique(fst$pop1)

### sort out and remove dups
fst$poppair <- paste(fst$pop1,fst$pop2,sep = "_vs_")
test <- c(fst$poppair,paste(fst$pop2,fst$pop1,sep = "_vs_"))
dups <- names(table(test)[table(test)>1])
rem.lsit <- c()
for (i in 1:nrow(fst)) {
  fst[i,1:2] <- sort(c(fst[i,1],fst[i,2]))
  fst$poppair[i] <- paste(fst$pop1[i],fst$pop2[i],sep = "_vs_")
  if(fst$pop1[i]!=fst$pop2[i] & any(dups == paste(fst$pop2[i],fst$pop1[i],sep = "_vs_"))){
    if(!any(paste(fst$pop1[rem.lsit],fst$pop2[rem.lsit],sep = "_vs_") == paste(fst$pop1[i],fst$pop2[i],sep = "_vs_"))) {
      rem.lsit <- c(rem.lsit,i)
    }
  }
}

fst <- fst[-rem.lsit,]
colnames(fst) <- c("pop1","pop2","Fst","poppair")


#### add geodistance
for (i in 1:nrow(fst)) {
  fst$dis2[i] <- geodist(info[which(info$pop == (fst$pop1[i]))[1],6:7],
                         info[which(info$pop == (fst$pop2[i]))[1],6:7],measure = "geodesic")/1000
}


####################################### cnv distance, weight manhatan distance
cnv.dep.final.noscal.popmean <- apply(cnv.dep.final,1,function(x){
  tapply(x, info$pop, mean)
})
probe.mean <- rowMeans(cnv.dep.final)
cnv.dep.final.noscal.popmean <- cnv.dep.final.noscal.popmean[-1,]

fst$cnv_dis2 <- apply(fst,1,FUN=function(x){
  probe_dis <- apply(rbind(cnv.dep.final.noscal.popmean[c(x[1],x[2]),],probe.mean),2, 
                     FUN = function(x) {
                       abs(x[1]-x[2])/x[3]
                     })
  return(sum(probe_dis))
})

par(mfrow = c(1,3))
plot(fst$cnv_dis2,fst$Fst)
plot(fst$cnv_dis,fst$Fst)
plot(fst$cnv_dis2,fst$cnv_dis)

############################################# env distance 
#### extract env data
env <- data.frame(indi=info$indi,pop=info$pop,lat=info$Lat,long=info$Long)

for (i in 1:19){
  bio_raster<- raster(paste("CHELSA_bio",i,"_1981-2010_V.2.1.tif",sep = ""))
  env[,paste("bio",i,sep = "")] <- extract(bio_raster,sw.env.indi[,c(4,3)])
}
saveRDS(env,"env.rds")

#### distance calulation
env <- apply(env[,-(1:2)],2,function(x){tapply(as.numeric(x), info$pop[match(env.cline$indi,info$indi)], mean)})
env.dist <- as.matrix(dist(scale(env[,3:21])))
colnames(env.dist) <- sort(rownames(env))
rownames(env.dist) <- sort(rownames(env))

for (i in 1:nrow(fst)) {
  fst$envdist[i] <- env.dist[fst$pop1[i],fst$pop2[i]]
}

#### mentel test for IBD
fst_mat <- acast(fst, pop1 ~ pop2, value.var = "cnv_dis2")
fst_mat[lower.tri(fst_mat)] <- t(fst_mat)[lower.tri(fst_mat)]  # make symmetric
fst_dist <- as.dist(fst_mat)

dist_mat <- acast(fst, pop1 ~ pop2, value.var = "dis2")
dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
geo_dist <- as.dist(dist_mat)

mantel(fst_dist, geo_dist, method = "pearson", permutations = 9999)


### pairwise correlation
cor.test(fst$Fst,fst$cnv_dis2,method = "spearman")
cor.test(fst$Fst,fst$dis2,method = "spearman")
cor.test(fst$Fst,fst$envdist,method = "spearman")
cor.test(fst$dis2,fst$cnv_dis,method = "spearman")
cor.test(fst$envdist,fst$cnv_dis2,method = "spearman")
cor.test(fst$envdist,fst$dis2,method = "spearman")

cor.test(residuals(lm(fst$Fst~fst$dis2)),fst$envdist,method = "spearman")
cor.test(residuals(lm(fst$cnv_dis~fst$dis2)),fst$envdist,method = "spearman")
cor.test(residuals(lm(fst$envdist~fst$dis2)),fst$cnv_dis,method = "spearman")


############# plot
#  Fst ~ dis 
cor.test(fst$Fst,fst$dis2,method = "spearman")
ggplot(fst,aes(dis2,Fst))+ geom_point(color = brewer.pal(6,"Greens")[4]) + 
  xlab("Geodesic distance (km)")+
  ylab("Fst (Single-copy SNP)")+
  geom_smooth(method=lm , color = brewer.pal(6,"GnBu")[6], linetype = "dashed", se=TRUE) +
  theme_classic()+
  annotate(
    "text", 
    x = 4000, 
    y = 0.1, 
    label = "Spearman's rho = 0.90
    p < 0.001",
    size = 5,
    color = "black"
  ) + ylim(0,0.36)+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1),
        plot.margin = margin(l = 14,r = 4))


# cnv_dis ~ dis
cor.test(fst$cnv_dis2,fst$dis2,method = "spearman")
ggplot(fst,aes(dis2,cnv_dis2))+ geom_point(color = brewer.pal(6,"Greens")[4]) + 
  xlab("Geodesic distance (km)")+
  ylab("CNV distance")+
  geom_smooth(method=lm , color = brewer.pal(6,"GnBu")[6], linetype = "dashed", se=TRUE) +
  theme_classic()+
  annotate(
    "text", 
    x = 4000, 
    y = 500, 
    label = "Spearman's rho = 0.67
    p < 0.001",
    size = 5,
    color = "black"
  ) + ylim(250,1000)+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))



# cnv_dis ~ fst
cor.test(fst$cnv_dis2,fst$Fst,method = "spearman")
ggplot(fst,aes(Fst,cnv_dis2))+ geom_point(color = brewer.pal(6,"Greens")[4]) + xlab("Fst (Single-copy SNP)")+
  ylab("CNV distance")+
  geom_smooth(method=lm , color = brewer.pal(6,"GnBu")[6], linetype = "dashed", se=TRUE,) +
  theme_classic()+
  annotate(
    "text", 
    x = 0.25, 
    y = 500, 
    label = "Spearman's rho = 0.70
    p < 0.001",
    size = 5,
    color = "black"
  )+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))


### pc correlation between SNP and gCNV
cor.test(evec$PC1,cnv.dep.final.rotation$PC1,method = "spearman")
ggplot(data.frame(SNP = evec$PC1, CNV = cnv.dep.final.rotation$PC1),aes(SNP,CNV))+ 
  geom_point(color = brewer.pal(6,"Greens")[4]) + xlab("PC1 (Single-copy SNP)")+
  ylab("PC1 (CNV)")+
  geom_smooth(method=lm , color = brewer.pal(6,"GnBu")[6], linetype = "dashed") +
  theme_classic()+
  annotate(
    "text", 
    x = 0.06, 
    y = 8, 
    label = "Spearman's rho = 0.95
    p < 0.001",
    size = 5,
    color = "black"
  )+ 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1),
        plot.margin = margin(l = 14,r=6))

