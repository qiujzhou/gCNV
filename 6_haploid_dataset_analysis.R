library(MASS)
library(rCNV)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)
library(ggpubr)
library(RColorBrewer)
library(ggsignif)
num_cores = 20
cl <- makeCluster(num_cores)

### basic statistics
mega.info <- read.table("hap.info.txt",h=T)  #### sample info table: ID, pop, lat, long..


#### sample heterozygosity distribution
sw.het <- read.table("filepath",h=T)    ### heterozygosity output from vcftools
cline.het <- read.table("filepath",h=T)
mega.het <- read.table("filepath",h=T)

sw.het <- separate(sw.het, OBS.HOM1.HET.HOM2., into = c("H1", "HET", "H2"), sep = "/",convert = T)
cline.het <- separate(cline.het, OBS.HOM1.HET.HOM2., into = c("H1", "HET", "H2"), sep = "/",convert = T)
mega.het <- separate(mega.het, OBS.HOM1.HET.HOM2., into = c("H1", "HET", "H2"), sep = "/",convert = T)

sw.het$heter <- sw.het$HET/(sw.het$H1+sw.het$HET+sw.het$H2) 
cline.het$heter <- cline.het$HET/(cline.het$H1+cline.het$HET+cline.het$H2) 
mega.het$heter <- mega.het$HET/(mega.het$H1+mega.het$HET+mega.het$H2) 

mean(sw.het$het)
mean(cline.het$het)
mean(mega.het$het)
sd(sw.het$het)
sd(cline.het$het)
sd(mega.het$het)


boxplot(c(sw.het$het,cline.het$het,mega.het$het) ~ c(rep("Swedish cline",nrow(sw.het)),
                                                     rep("P.abies-P.obovata",nrow(cline.het)),
                                                     rep("Megagametophyte",nrow(mega.het))))
het <- data.frame(het = c(sw.het$het,cline.het$het,mega.het$het),dataset = c(rep("Swedish cline",nrow(sw.het)),
                                                                             rep("P.abies-P.obovata",nrow(cline.het)),
                                                                             rep("Megagametophyte",nrow(mega.het))))

ggplot(het, aes(x = dataset, y = het)) +
  #geom_violin(width = 2,trim = FALSE,col="darkgrey") + ylim(0,0.4) + 
  theme_minimal() + 
  geom_boxplot(width = 0.5, fill = "white",col="darkgray", outlier.shape = NA) +
  theme(legend.position = "none") +
  theme_classic() + 
  labs(title = "Observed heterozygosity", x = "", y = "Observed heterozygosity")+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


##### hap-based validation
mega.gt <- readRDS("mega.gt.rds") ### hap genotype table

#### heterozygosity per probe
mega.site.het <- parApply(mega.gt[,-(1:4)],1,FUN = function(x){
  sum(x == "0/1" | x == "1/0")/length(x[x!="./."])
},cl=cl)

mega.snp.info <- data.frame(mega.gt[,1:4],het = mega.site.het)

#### add probe and gene
mega.snp.info$probe <- parApply(mega.snp.info,1,function(x, probe ){
  if (any(probe$V4 < as.numeric( x[2]) & probe$V5 > as.numeric(x[2]) & probe$V1 == as.character(x[1])) ) {
    return(probe$probe[probe$V4 < as.numeric( x[2]) & probe$V5 > as.numeric(x[2]) & probe$V1 == as.character(x[1])][1])
  } else{
    return(FALSE)
  }
},cl = cl,probe=gff)

mega.snp.info$gene <- parApply(mega.snp.info,1,function(x, probe ){
  if (any(probe$V4 < as.numeric( x[2]) & probe$V5 > as.numeric(x[2]) & probe$V1 == as.character(x[1])) ) {
    return(probe$gene[probe$V4 < as.numeric( x[2]) & probe$V5 > as.numeric(x[2]) & probe$V1 == as.character(x[1])][1])
  } else{
    return(FALSE)
  }
},cl = cl,probe=gff)
mega.snp.info$gene <- probe2gene$gids[match(mega.snp.info$gene,probe2gene$gene)]

#### mean het per gene
probe.het <- tapply(mega.snp.info$het, mega.snp.info$gene, mean,na.rm=T)
probe.het <- probe.het[-1]

probe.info <- data.frame(gene = unique(probe2gene$gids),het = probe.het[unique(probe2gene$gids)])
#probe.info$dep <- tapply(mega.snp.info$dep, mega.snp.info$probe, mean,na.rm=T)[-1]


### het difference between CNV and non-CNV
probe.info$sw_stat <- "Unknown"
probe.info$sw_stat[probe.info$gene %in% sw.cnv.gene] <- "gCNV"
probe.info$sw_stat[probe.info$gene %in% sw.nonCNV.gene] <- "Single copy gene"
probe.info$er_stat <- "Unknown"
probe.info$er_stat[probe.info$gene %in% cline.cnv.gene] <- "gCNV"
probe.info$er_stat[probe.info$gene %in% cline.nonCNV.gene] <- "Single copy gene"
probe.info$sw_dep <- tapply(sw.probe.stat$cov,sw.probe.stat$gene,mean)[probe.info$gene]
probe.info$er_dep <- tapply(cline.probe.stat$cov,cline.probe.stat$gene,mean)[probe.info$gene]
probe.info$sw_cv <- tapply(sw.probe.stat$cv,sw.probe.stat$gene,mean)[probe.info$gene]
probe.info$er_cv <- tapply(cline.probe.stat$cv,cline.probe.stat$gene,mean)[probe.info$gene]

probe.info2 <- probe.info %>% 
  gather(key = "dataset", value = "Group", sw_stat, er_stat)
probe.info2$dataset[probe.info2$dataset == "sw_stat"] <- "Swedish cline"
probe.info2$dataset[probe.info2$dataset == "er_stat"] <- "P. abies-P. obovata"
probe.info2$dep[probe.info2$dataset == "Swedish cline"] <- probe.info2$sw_dep[probe.info2$dataset == "Swedish cline"]
probe.info2$dep[probe.info2$dataset == "P. abies-P. obovata"] <- probe.info2$er_dep[probe.info2$dataset == "P. abies-P. obovata"]
probe.info2$CV[probe.info2$dataset == "Swedish cline"] <- probe.info2$sw_cv[probe.info2$dataset == "Swedish cline"]
probe.info2$CV[probe.info2$dataset == "P. abies-P. obovata"] <- probe.info2$er_cv[probe.info2$dataset == "P. abies-P. obovata"]
probe.info2$Group <- as.factor(probe.info2$Group)
probe.info2$dataset <- as.factor(probe.info2$dataset)

####plot
ggplot(probe.info2[!probe.info2$Group %in% "Unknown",], aes(x = dataset, y = het, col=Group)) +
  theme_minimal() + ylim(0,0.2) + 
  geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA,alpha = 0.5) +
  theme(legend.position = "none") +
  theme_classic() + 
  labs(title = "Observed heterozygosity", x = "", y = "Observed heterozygosity")+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))+
  annotate("text", x = 1, y = 0.15, label = "***",size = 5)+
  annotate("text", x = 2, y = 0.15, label = "***",size = 5)+
  #theme(axis.text.x = element_text(angle = 30, hjust = 1))+ 
  scale_color_manual(values= c(brewer.pal(8,"Set2"))[c(1,2)])

length(unique(probe.info$gene[probe.info$sw_stat %in% "gCNV"]))
length(unique(probe.info$gene[probe.info$sw_stat %in% "Single copy gene"]))
length(unique(probe.info$gene[probe.info$er_stat %in% "gCNV"]))
length(unique(probe.info$gene[probe.info$er_stat %in% "Single copy gene"]))

#### statistical test
t.test(probe.info$het[!probe.info$sw_stat %in% "Unknown"]~
         probe.info$sw_stat[!probe.info$sw_stat %in% "Unknown"])
t.test(probe.info$het[!probe.info$er_stat %in% "Unknown"]~
         probe.info$er_stat[!probe.info$er_stat %in% "Unknown"])
sd(probe.info$het[probe.info$sw_stat %in% "gCNV"],na.rm = T)
sd(probe.info$het[probe.info$sw_stat %in% "Single copy gene"],na.rm = T)
sd(probe.info$het[probe.info$er_stat %in% "gCNV"],na.rm = T)
sd(probe.info$het[probe.info$er_stat %in% "Single copy gene"],na.rm = T)
