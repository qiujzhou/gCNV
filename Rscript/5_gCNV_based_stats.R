library(rCNV)
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(parallel)
library(ggpubr)


################################################### check probe depth distribution and CV
probe.stat <- data.frame(probe = rownames(depth_final))
probe.stat$cov <- rowMeans(depth_final)
probe.stat$cv <- apply(depth_final,1,FUN=function(x){
  sd(x,na.rm = T)/mean(x,na.rm = T)
})
probe.stat$stat <- NA
probe.stat$stat[probe.stat$probe %in% cnv.probe.final] <- "CNV"
probe.stat$stat[probe.stat$probe %in% nonCNV.probe] <- "Single copy"
probe.stat$gene <- probe2gene$gids[match(probe.stat$probe,probe2gene$probe)]

sd(probe.stat$cov[probe.stat$stat == "CNV"],na.rm = T)
sd(probe.stat$cov[probe.stat$stat == "Single copy"],na.rm = T)
sd(probe.stat$cv[probe.stat$stat == "CNV"],na.rm = T)
sd(probe.stat$cv[probe.stat$stat == "Single copy"],na.rm = T)
t.test(probe.stat$cov[probe.stat$stat == "CNV"],
       probe.stat$cov[probe.stat$stat == "Single copy"])
t.test(probe.stat$cv[probe.stat$stat == "CNV"],
       probe.stat$cv[probe.stat$stat == "Single copy"])


#################################################### pi vs allele size variation
alv <- function(dep.tab,pop,cl){
  alv_probe_mat <- parallel::parApply(dep.tab,1,function(x,pop){
    alv_probe <- tapply(x, pop, function(x){
      var(x)
    })
    return(alv_probe)
  },cl=cl,pop=pop)
  return(unlist(rowSums(alv_probe_mat)))
}

alv.test1 <- alv(cnv.dep.final,pop = info$pop, cl = cl)

pop.pi <- read.table("summary.pi") #### SNP-based pi from vcftools
#ref.gff <- read.table("/home/zhouqj/cline_proj/ref/Pabies.only.gene.gff3")
colnames(pop.pi) <- c("Pop","Nsite","pi_sum")
pop.pi$pi <- pop.pi$pi_sum/pop.pi$Nsite


cor.test(pop.pi$pi,pop.pi$weighta.alv,method = "spearman")
p.pi.alv <- ggplot(pop.pi,aes(pi,weighta.alv))+ geom_point(color = brewer.pal(6,"Greens")[4]) + 
  xlab("Pi (Single-copy SNP)")+
  ylab("Allele size variance")+
  theme_classic()+
  annotate(
    "text", 
    x = 0.14, 
    y = 3500, 
    label = "Spearman's rho = -0.24
    p = 0.09",
    size = 5,
    color = "black"
  )+
  geom_smooth(method=lm , color = brewer.pal(6,"GnBu")[6], linetype = "dashed", se=TRUE) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1),
        plot.margin = margin(l = 14,r = 4))
p.pi.alv



