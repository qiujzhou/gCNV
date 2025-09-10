library(rCNV)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(parallel)
library(ggpubr)
library(mclust)



### pca on nonCNV probes
noncnv.pca <- prcomp(depth_final[,colnames(depth_final) %in% nonCNV.probe],
                                scale. = T,center = T,retx = T)
noncnv.rotation <- cbind(noncnv.pca$x[,1:10],
                                    info[match(rownames(depth_final),info$indi),])
noncnv.pca$sdev[1:2]^2/sum(noncnv.pca$sdev^2) #### variance explianed by PC1,2
rm(noncnv.pca)
#### color by genetic Latitude
p3<- ggplot(noncnv.rotation) + geom_point(aes(PC1,PC2,col=Long))+
  xlab("PC1 (3.3%)") + ylab("PC2 (1.9%)") + 
  theme_classic()+
  labs(color = "Longitude")+
  ggtitle("Single-copy gene probes") +
  scale_color_continuous(low = "purple", high = "orange") +
  theme(
    panel.grid.major = element_blank(),    # Remove major gridlines
    panel.grid.minor = element_blank() 
  ) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))
p3

#### color by genetic cluster
p3.cluster <- ggplot(noncnv.rotation) + geom_point(aes(PC1,PC2,col=cluster))+
  xlab("PC1 (3.3%)") + ylab("PC2 (1.9%)") + 
  theme_classic()+
  scale_color_manual(values= c(brewer.pal(9,"Set1"))[c(2:6,9,7,8,1)]) + 
  labs(color="Genetic cluster") +
  ggtitle("Single-copy gene probes") +
  theme(
    panel.grid.major = element_blank(),    # Remove major gridlines
    panel.grid.minor = element_blank() 
  ) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))
p3.cluster

### pca on CNV probes
cnv.dep.final <- depth_final[,intersect(cnv.probe.final,colnames(depth_final))]
cnv.dep.final.pca <- prcomp(cnv.dep.final,scale. = F,center = F,retx = T)
cnv.dep.final.rotation <- cbind(cnv.dep.final.pca$x[,1:10],info[match(rownames(cnv.dep.final),info$indi),])
cnv.dep.final.pca$sdev[1:2]^2/sum(cnv.dep.final.pca$sdev^2)


##colored by cluster
p4.cluster <- ggplot(cnv.dep.final.rotation) + geom_point(aes(PC1,PC2,col=cluster))+
  xlab("PC1 (5.0%)") + ylab("PC2 (2.2%)") + 
  theme_classic()+
  scale_color_manual(values= c(brewer.pal(9,"Set1"))[c(2:6,9,7,8,1)]) + 
  labs(color="Genetic cluster") +
  ggtitle("CNV probes") +
  theme(
    panel.grid.major = element_blank(),    # Remove major gridlines
    panel.grid.minor = element_blank() 
  ) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))
p4.cluster

#colored by Longitude
p4 <- ggplot(cnv.dep.final.rotation) + geom_point(aes(PC1,PC2,col=Long))+
  xlab("PC1 (4.6%)") + ylab("PC2 (2.2%)") + 
  theme_classic()+
  #  scale_color_manual(values= c(brewer.pal(9,"Set1"))[c(2:6,9,7,8,1)]) + 
  labs(color="Longitude") +
  ggtitle("CNV probes") +
  scale_color_continuous(low = "purple", high = "orange") +
  theme(
    panel.grid.major = element_blank(),    # Remove major gridlines
    panel.grid.minor = element_blank() 
  ) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))
p4


