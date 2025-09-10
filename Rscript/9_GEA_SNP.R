library(raster)
library(sp)
require(rworldmap)
require(car)
require(raster)
require(ade4)
require(sp)
#require(readxl)
require(corrplot)
require(lme4)
require(lmtest)
require(lmerTest)
library(fitdistrplus)
require(geosphere)
require(FactoMineR)
require(factoextra)
require(gridExtra)
library(gridGraphics)
#install.packages("fitdistrplus")
library("fitdistrplus")


################################################################### 
###### prepare input for bayenv2
#saveRDS(freq.true.snp,"freq.true.snp.rds")
all_freq <- fread("frq.strat",h=T)  #### strat output from PLINK
all_freq$MaAC <- all_freq$NCHROBS - all_freq$MAC


##### for baypass
all_freq1 <- all_freq[,c(3,7,9)] %>%  group_by(CLST) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = CLST, values_from = c(MAC,MaAC)) %>%
  dplyr::select(-row)

snp.pos <- fread("snp.pos") ### all SNP
colnames(snp.pos) <- c("CHR","POS")
snp.pos$index <- paste(snp.pos$CHR,snp.pos$POS,sep = "_")
cnv.pos <- read.table("snp.cnv.pos")  ### SNP within CNV
colnames(cnv.pos) <- c("CHR","POS")
cnv.pos$index <- paste(cnv.pos$CHR,cnv.pos$POS,sep = "_")

SNP_filter <- snp.pos$index %in% cnv.pos$index  ### filter
all_freq1 <- all_freq1[SNP_filter,]  #### for SNP based glm

##### for baypass
order_allele <- c(1:110)  ### 55 populations
order_allele[(1:55)*2-1] <- 1:55
order_allele[(1:55)*2] <- (1:55)+55
all_freq1 <- as.data.frame(all_freq1)[,c(order_allele)]
write.table(all_freq1,"freq_for_baypass_CNVsnp.txt",quote = F,col.names = F, row.names = F)

rm(all_freq)

##### for bayenv2
all_freq <- parApply(all_freq1, 1, FUN = function(x){
  temp <- matrix(nrow = 2,ncol = 55)
  temp[1,] <- x[1:55*2-1]
  temp[2,] <- x[1:55*2]
  return(data.frame(temp))
},cl=cl)
all_freq <- do.call(rbind,all_freq)
write.table(all_freq,"all_freq_bayenv2.txt", col.names = F,row.names = F,quote = F,sep = "\t")
write.table(t(round(scale(env.pop[,c(1:19,25:27)]),2)),"env.txt",col.names = F,row.names = F,quote = F,sep = "\t")


###################### glm on SNPs 
env.pop <- apply(env[,c(5:31)], 2, function(x){
  tapply(x, env$pop, mean)
})

env.pop <- as.data.frame(env.pop)

gea.snp.p <- pbapply(all_freq1, 1, FUN=function(x,env){
  stat <- c()
  for (i in c(1:19,25:27)) {
    env2 <- data.frame(mac=x[1:55],mic = x[56:110],bio = env[,i], 
                       PC1_SNP = env$PC1_SNP,
                       PC2_SNP = env$PC2_SNP,
                       PC3_SNP = env$PC3_SNP,
                       PC4_SNP = env$PC4_SNP,
                       PC5_SNP = env$PC5_SNP
    )
    test <- summary(glm(cbind(mac,mic) ~ scale(bio)+PC1_SNP+PC2_SNP+PC3_SNP+PC4_SNP+PC5_SNP, 
                        family = binomial,
                        data = env2 ))
    stat <- c(stat,test$coefficients[2,4])
  }
  names(stat) <- colnames(env)[c(1:19,25:27)]
  return(stat)
},cl=cl,env=as.data.frame(env.pop))


cnv.pos <- fread("snp.cnv.pos")
colnames(cnv.pos) <- c("CHR","POS")
cnv.pos$probe <- parApply(cnv.pos,1,function(x, probe ){
  if (any(probe$V4 < as.numeric( x[2]) & probe$V5 > as.numeric(x[2]) & probe$V1 == as.character(x[1])) ) {
    return(probe$probe[probe$V4 < as.numeric( x[2]) & probe$V5 > as.numeric(x[2]) & probe$V1 == as.character(x[1])][1])
  } else{
    return(FALSE)
  }
},cl = cl,probe=gff)

##### adjusted p value
gea.snp.p.adj <- gea.snp.p
for (i in 1:22) {
  gea.snp.p.adj[i,] <- p.adjust(gea.snp.p[i,])
}

##### keep only SNP with min p for each probe 
gea.snp.p.min <- do.call(cbind,lapply(cnv.probe.filter, function(x){
  apply(gea.snp.p[,cnv.pos$probe %in% x], 1, min)
}))

gea.snp.p.adj.min <- do.call(cbind,lapply(cnv.probe.filter, function(x){
  apply(gea.snp.p.adj[,cnv.pos$probe %in% x], 1, min)
}))

for (i in 1:22) {
  print( sum(gea.snp.p.adj.min[i,] < 0.2))
}

### plot
plot_merge <- list()
for (i in 1:22){
  dat <- data.frame(x = -log10(gea.cnv.p[i,]),
                    y = -log10(gea.snp.p.min[i,]),
                    category = "Not significant")
  dat$category[p.adjust(gea.cnv.p[i,]) < 0.2] <- "Significant in gCNV based GEA"
  dat$category[gea.snp.p.adj.min[i,] < 0.2] <- "Significant in SNP based GEA"
  dat$category[gea.snp.p.adj.min[i,] < 0.2 & p.adjust(gea.cnv.p[i,]) < 0.2] <- 
    "Significant in gCNV & SNP based GEA"
  dat$category <- factor(dat$category, levels = c(
    "Not significant",
    "Significant in gCNV based GEA",
    "Significant in SNP based GEA",
    "Significant in gCNV & SNP based GEA"
  ))
  
  dat <-  dat[is.finite(dat$y),]
  plot_merge[[i]] <- ggplot(data = dat,aes(x,y)) + 
    geom_point(aes(color = category)) + 
    labs(
      #y = expression(-log[10](italic(p))~ (SNP_GLM)),
      #x = expression(-log[10](italic(p))~ (gCNV_GLM)), 
      x = NULL, y =NULL,
      title = rownames(gea.cnv.p)[i]) +
    geom_smooth(method=lm , color = "#117733", linetype = "dashed", se=TRUE,) +
    theme_classic()+
    annotate(
      "text", 
      x = max(dat$x)/2, 
      y = max(dat$y)*0.6, 
      label = paste("rho =",round(cor(-log10(gea.cnv.p[i,]),-log10(gea.snp.p.min[i,]),method = "spearman"),2)),
      size = 3,
      color = "black"
    )+theme(plot.title = element_text(hjust = 0.5, vjust = 1)) + 
    scale_color_manual(drop = T,
                       values = c(
                         "Not significant" = "#B0B0B0",       
                         "Significant in gCNV based GEA" = "#FF7F00", 
                         "Significant in SNP based GEA" = "#1F78B4",
                         "Significant in gCNV & SNP based GEA" = "#984EA3"
                       )) + scale_y_continuous(breaks = seq(0, 8, 2),
                                               limits = c(0,max(dat$y)+0.5)) 
}


for (i in 1:22){
  print(length(intersect(colnames(gea.cnv.p)[p.adjust(gea.cnv.p[i,]) < 0.2],
                         colnames(gea.cnv.p)[gea.snp.p.adj.min[i,]< 0.2] )))
}

gea.cnv.probe.bio <- c()
for (i in 1:19){
  gea.cnv.probe.bio <- c(gea.cnv.probe.bio,
                         colnames(gea.cnv.p)[p.adjust(gea.cnv.p[i,]) < 0.2])
}
length(unique(probe2gene$gids[probe2gene$probe %in% gea.cnv.probe.bio]))


gea.cnv.probe.pc <- c()
for (i in 20:22){
  gea.cnv.probe.pc <- c(gea.cnv.probe.pc,
                        colnames(gea.cnv.p)[p.adjust(gea.cnv.p[i,]) < 0.2])
}
length(unique(probe2gene$gids[probe2gene$probe %in% gea.cnv.probe.pc]))

length(unique(probe2gene$gids[probe2gene$probe %in% 
                                intersect(gea.cnv.probe.pc,gea.cnv.probe.bio)]))



#### results from bayenv2
for (i in 1:10) {
  factor19_bf <- read.table(paste("bayenv2/new/batch",i,"/bf_environ.env.txt",sep = ""))
  factor19_bf <- factor19_bf[order(factor19_bf$V1),]
  write.table(factor19_bf,paste("bayenv2/new/batch",i,"/bf2_environ.env.txt",sep = ""),
              col.names = F,row.names = F,quote = F)
}


factor19_bf <- fread("all.bf.txt",h=F)
factor19_bf <- cbind(cnv.pos,factor19_bf[,-1])

##### keep only SNP with the maximum bayes factor
bay.snp.bf.max <- do.call(cbind,lapply(colnames(gea.cnv.p), function(x){
  apply(factor19_bf[factor19_bf$probe %in% x,-(1:4)], 2, function(x){
    max(abs(as.numeric(x)))
  })
}))
colnames(bay.snp.bf.max) <- colnames(gea.cnv.p)


plot_merge <- list()
for (i in 1:22){
  dat <- data.frame(x = -log10(gea.cnv.p[i,]),
                    y = bay.snp.bf.max[i*3-2,],
                    category = "Not significant")
  dat$category[p.adjust(gea.cnv.p[i,]) < 0.2] <- "Significant in gCNV based GEA"
  dat$category[dat$y > 20] <- "Significant in SNP based GEA"
  dat$category[dat$y > 20 & p.adjust(gea.cnv.p[i,]) < 0.2] <- 
    "Significant in gCNV & SNP based GEA"
  dat$category <- factor(dat$category, levels = c(
    "Not significant",
    "Significant in gCNV based GEA",
    "Significant in SNP based GEA",
    "Significant in gCNV & SNP based GEA"
  ))
  
  plot_merge[[i]] <- ggplot(dat,aes(x,y)) + 
    geom_point(aes(color = category)) +
    labs(x = NULL, y = NULL, title = rownames(sw.gea.cnv.p)[i]) +
    geom_smooth(method=lm , color = "#117733", linetype = "dashed", se=TRUE,) +
    theme_classic()+
    annotate(
      "text", 
      x = max(dat$x)/2, 
      y = max(dat$y)*0.6, 
      label = paste("rho =",round(cor(-log10(gea.cnv.p[i,]),bay.snp.bf.max[i*3-2,],method = "spearman"),2)),
      size = 3,
      color = "black"
    )+theme(plot.title = element_text(hjust = 0.5, vjust = 1)) +
    scale_color_manual(drop = FALSE,
                       values = c(
                         "Not significant" = "#B0B0B0",       
                         "Significant in gCNV based GEA" = "#FF7F00", 
                         "Significant in SNP based GEA" = "#1F78B4",
                         "Significant in gCNV & SNP based GEA" = "#984EA3"
                       )) 
}


for (i in 1:22){
  #print(sum(bay.snp.bf.max[i*3-2,]>20))
  print(length(intersect(colnames(bay.snp.bf.max)[bay.snp.bf.max[i*3-2,]>20],
                         colnames(gea.cnv.p)[p.adjust(gea.cnv.p[i,]) < 0.2])))
}
