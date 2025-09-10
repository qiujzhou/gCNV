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



env <- readRDS("env.rds")
#### add PC from SNP data
env$PC1_SNP <- info$pc1[match(env$indi,info$indi)]
env$PC2_SNP <- info$pc2[match(env$indi,info$indi)]
env$PC3_SNP <- info$pc3[match(env$indi,info$indi)]
env$PC4_SNP <- info$pc4[match(env$indi,info$indi)]
env$PC5_SNP <- info$pc5[match(env$indi,info$indi)]

colnames(env)
### environmental variable exploration
res.env<-PCA(env[,c(3:23)],quanti.sup = c(1,2),graph = T)
eigenvalues <- sw.res.env$eig
eigenvalues
plot.var.env<-fviz_pca_var(sw.res.env, axes = c(1,2),col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = TRUE , title= "Environmnetal data")
plot(res.env$ind$coord[,2]~res.env$ind$coord[,1],pch=19)
plot.var.env
eigenvalues <- res.env$eig
eigenvalues
plot.var.env<-fviz_pca_var(res.env, axes = c(1,2),col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = TRUE , title= "Environmnetal data") + 
  theme(
    plot.margin = margin(15, 25, 15, 15), # top, right, bottom, left
    plot.title = element_text(
      face = "bold",        # bold
      hjust = 0.5,          # center horizontally
      size = 14             # optional: increase size
    )
  )
plot.var.env
plot(res.env$ind$coord[,2]~res.env$ind$coord[,1],pch=19)
var.env <- get_pca_var(res.env)
plot.cor.env<-corrplot(var.env$contrib, is.corr=FALSE,tl.col = "black",tl.cex = 0.6,cl.pos = 'n',)
plot.cor.env
plot.scre.env<-fviz_screeplot(res.env, ncp=6,main="Screeplot - Environmental data",
                              xlab="Principal component (PC)")


#### association between top PCs from ENV, SNP
head(info)
info$ePC1<-res.env$ind$coord[,1]
info$ePC2<-res.env$ind$coord[,2]
info$ePC3<-res.env$ind$coord[,3]
env$ePC1 <- res.env$ind$coord[,1]
env$ePC2 <- res.env$ind$coord[,2]
env$ePC3 <- res.env$ind$coord[,3]
env$gPC1 <- cnv.dep.final.rotation$PC1
env$gPC2 <- cnv.dep.final.rotation$PC2
env$gPC3 <- cnv.dep.final.rotation$PC3

outGE<-matrix(NA,ncol = 3, nrow = 3)
rownames(outGE)<-c("gDim.1","gDim.2","gDim.3")
colnames(outGE)<-c("eDim.1","eDim.2","eDim.3")

poutGE<-matrix(NA,ncol = 3, nrow = 3)
rownames(poutGE)<-c("gDim.1","gDim.2","gDim.3")
colnames(poutGE)<-c("eDim.1","eDim.2","eDim.3")
for(i in c(1:3)){
  for(j in c(1:3)){
    outGE[i,j]<-abs(cor(env[,28+i],env[,31+j],method = "spearman"))
    
    poutGE[i,j]<-cor.test(env[,28+i],env[,31+j],method = "spearman")$p.value
    
  }
}


plot.var.env
corrplot(outGE,
         main = "Genotype - Environment",
         tl.col = "black", tl.srt = 45, p.mat = poutGE,
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 2,
         insig = 'label_sig', pch.col = 'grey5',
         cl.align.text = 'l', mar = c(2,0,2,2),
         tl.cex = 1.1, cex.main = 1,
         cl.length = 5, is.corr = FALSE)

grid.echo()
corr_grob <- grid.grab()

ggarrange(plot.var.env, corr_grob, ncol = 2,
                   widths = c(1.1,1),labels = c("A","B"))



############################# GEA analysis
###### PCadapt, explore number of PCs used for controlling pop structure
library(pcadapt)
library(VennDiagram)
bed_SNPs <- read.pcadapt("path2bed",
                                   type = "bed")
sw.x <- pcadapt(input = bed_SNPs, K = 20) 
plot(sw.x, option = "screeplot")  ### best k = 3


gea.cnv.p <- parApply(depth_final[cnv.probe.final,], 1, FUN=function(x,env){
  stat <- c()
  for (i in c(5:23,29:31)) {  ### 19 bioclimatic variables, and top 3 ePCs
    env2 <- data.frame(dep=x,pop = env[,2], bio = env[,i],env[,c(24:28)])
    test <- summary(glm(dep ~ scale(bio) + PC1_SNP + PC2_SNP + PC3_SNP + PC4_SNP + PC5_SNP, 
                        family =  "gaussian", data = env2 ))
    #test <- summary(lm(dep ~ scale(bio) , data = env2 ))
    #test <- cor.test(env2$dep, scale(env2$bio),data = env2,method="spearman")
    #stat <- c(stat,test$coefficients[2,4])
    stat <- c(stat,test$coefficients[2,4])
  }
  names(stat) <- colnames(env)[c(5:23,29:31)]
  return(stat)
},env=env,cl=cl)

### number of sig. probes & genes
for (i in 1:22) {
  sig.p <- colnames(gea.cnv.p)[p.adjust(gea.cnv.p[i,]) < 0.2]
  print(length(sig.p))
  print(length(unique(probe2gene$gids[probe2gene$probe %in% sig.p])))
}

### plot
plot_merge <- list()
for (i in 1:22) {
  dat <- data.frame(Observed = -log10(sort(gea.cnv.p[i,])),
                    Expected = -log10(ppoints(length(gea.cnv.p[i,]))),
                    Category = c("Not significant","Significant")[as.factor(p.adjust(sort(gea.cnv.p[i,])) < 0.2)])
  
  alpha <- 0.05
  n <- nrow(dat)
  dat <- dat %>%
    mutate(
      lower = -log10(qbeta(alpha / 2, seq_len(n), rev(seq_len(n)) + 1)),
      upper = -log10(qbeta(1 - alpha / 2, seq_len(n), rev(seq_len(n)) + 1)),
    )
  plot_merge[[i]] <- ggplot(dat, aes(Expected, Observed)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_point(aes(color = Category)) +
    scale_color_manual(values = c("grey40", "green3")) +
    theme_minimal() + ggtitle(rownames(gea.cnv.p)[i])+
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      panel.grid = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.25, 0.85),
      aspect.ratio = 1,    
      plot.title = element_text(
        face = "bold",        # bold
        hjust = 0.5,          # center horizontally
        size = 14             # optional: increase size
      )
    ) +
    labs(
      x = expression(Expected ~ -log[10](italic(p))),
      y = expression(Observed ~ -log[10](italic(p)))
    )
}

ggarrange(plotlist = plot_merge[1:19],ncol = 5,nrow = 4,
          common.legend = T,
          legend = "right")

