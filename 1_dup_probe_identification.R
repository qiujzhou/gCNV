library(MASS)
library(rCNV)
library(dplyr)
library(tidyr)
library(parallel)

##### parallelization and basic files
num_cores = 20
cl <- makeCluster(num_cores)
probe2gene <- readRDS("path2file") #### two columns table: 1st col: probe name,2nd:gene name
gff <- read.table("path2file",h=T) #### info table of probes: chr, start, end, probe name

##### processing raw vcf file
vcf <- readVCF("sample.vcf.gz")
vcf <- vcf$vcf
colnames(vcf)[1] <- "CHROM"

vcf.ad <- hetTgen(vcf)  #### allele depth table
vcf.gt <- hetTgen(vcf,"GT")  #### genotype table
vcf.tot <- parApply(vcf.ad[,-(1:4)],2,FUN=function(x){unlist(lapply(x,FUN=function(x){sum(as.numeric(unlist(strsplit(x,split = ","))))}))},cl=cl)
vcf.ad <- ad.correct(vcf.ad,gt.table = vcf.gt) #### ad correction

fis <- h.zygosity(vcf,cl=cl)
fis <- mean(fis$Fis,na.rm = T)


##### Run rCNV
alle.info <- allele.info(vcf.ad, Fis =fis, cl = cl,plot.allele.cov = F)
rcnv.dup <- cnv(alle.info, test=c("z.all","chi.all"), filter = "kmeans")
table(rcnv.dup$dup.stat)
table(rcnv.dup2$dup.stat)

##### SNP to probe
dup.info$probe <- parApply(dup.info,1,function(x, probe ){
  if (any(probe$start < as.numeric( x[2]) & probe$end > as.numeric(x[2]) & probe$chr == as.character(x[1])) ) {
    return(probe$probe[probe$start < as.numeric( x[2]) & probe$end > as.numeric(x[2]) & probe$chr == as.character(x[1])][1])
  } else{
    return(FALSE)
  }
},cl = cl,probe=gff)


true.snp <- rcnv.dup[rcnv.dup$dup.stat != "cnv",1:2]
true.cnv <- rcnv.dup[rcnv.dup$dup.stat == "cnv",1:2]

write.table(true.snp,"true.snp.pos",quote = F,row.names = F,col.names = F)
write.table(true.cnv,"true.cnv.pos",quote = F,row.names = F,col.names = F)


### probes enriched with putative cnv
snp_num <- tapply(dup.info$dup.stat,dup.info$probe,length)
enrich <- tapply(dup.info$dup.stat,dup.info$probe,FUN=function(x){
  y <- table(x)
  return(y["cnv"]/length(x))
  })
enrich[is.na(enrich)] <- 0


cnv_probe <- intersect(names(enrich[enrich>=0.3 ]), names( snp_num[snp_num * enrich>=2]))
noncnv_probe <- names(enrich[enrich<0.3 | snp_num * enrich <=1 ])[-1]
saveRDS(cnv_probe,"cnv.probe.rds")

### cnv_probe and noncnv_probe will be further refined using mclust

