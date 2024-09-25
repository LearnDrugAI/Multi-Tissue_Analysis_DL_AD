#### Make autoimmune disorder associated gene plots
## get subsets for each autoimmune disorder
library(data.table)
library(qvalue)

output_dir <- "Result_0.5"

# 检查目录是否存在，如果不存在则创建
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Directory created:", output_dir, "\n")
} else {
  cat("Directory already exists:", output_dir, "\n")
}


All_data_clean <- fread('All_data_clean.csv', header = TRUE, stringsAsFactors = FALSE, sep = ",")
datMHCmi_test <- All_data_clean
study_accessions <- c('GCST90041866', 'GCST90041879', 'GCST90042665', 'GCST90042678', 'GCST90042691')

qobjwo <- qvalue(datMHCmi_test$pvalue,fdr.level = 0.01)
datMHCmi_test$pva.qval <- qobjwo$qvalues
datMHCmi <- datMHCmi_test
datMHCmi$isSignificant <- qobjwo$significant

datMHCmi <- datMHCmi[pva.qval <= 0.5, ]
# datSigMHC <- datMHCmi[datMHCmi$isSignificant, ]
datSigMHC <- datMHCmi
#datSigMHC1=datSigMHC[abs(datSigMHC$zscore)>=mean(abs(datSigMHC$zscore)),]
datSigMHC1 <-datSigMHC


# to get associated genes with respect to each trait
TRAIT <- unique(datSigMHC$StudyAccession)
for (i in 1:length(TRAIT)){
  dattmp=subset(datSigMHC,datSigMHC$StudyAccession==TRAIT[i])
  dattmp=dattmp[,1:2]
  write.table(dattmp, 
              file = file.path(output_dir, paste0(TRAIT[i],'_AssoGeneMHC.txt')), 
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}

GCST90041866=subset(datSigMHC,datSigMHC$StudyAccession=='GCST90041866')
GCST90041879=subset(datSigMHC,datSigMHC$StudyAccession=='GCST90041879')
GCST90042665=subset(datSigMHC,datSigMHC$StudyAccession=='GCST90042665')
GCST90042678=subset(datSigMHC,datSigMHC$StudyAccession=='GCST90042678')
GCST90042691=subset(datSigMHC,datSigMHC$StudyAccession=='GCST90042691')

# ActDerm=subset(ActDerm,ActDerm$pred_perf_qval<=0.001)
# SLE=subset(SLE,SLE$pred_perf_qval<=0.001)
# CD=subset(CD,CD$pred_perf_qval<=0.001)
# UC=subset(UC,UC$pred_perf_qval<=0.001)
# RA=subset(RA,RA$pred_perf_qval<=0.001)

# GCST90041866=subset(GCST90041866,GCST90041866$pva.qval<=0.001)
# GCST90041879=subset(GCST90041879,GCST90041879$pva.qval<=0.001)
# GCST90042665=subset(GCST90042665,GCST90042665$pva.qval<=0.001)
# GCST90042678=subset(GCST90042678,GCST90042678$pva.qval<=0.001)
# GCST90042691=subset(GCST90042691,GCST90042691$pva.qval<=0.001)

# number of genes whose expression was associated with multiple traits
autoimmune=matrix(0,5,1)
autoimmune[1,1]=length(unique(GCST90041866$gene))
autoimmune[2,1]=length(unique(GCST90041879$gene))
autoimmune[3,1]=length(unique(GCST90042665$gene))
autoimmune[4,1]=length(unique(GCST90042678$gene))
autoimmune[5,1]=length(unique(GCST90042691$gene))
# make node size of autoimmune disorder
rownames(autoimmune)<-c('GCST90041866','GCST90041879','GCST90042665','GCST90042678','GCST90042691')
write.table(autoimmune, file = file.path(output_dir, "autoimmuneNodeSizeMHC.txt"), 
            row.names = TRUE, col.names = FALSE, quote = FALSE)
autoimmune <- read.table(file.path(output_dir, "autoimmuneNodeSizeMHC.txt"), 
                         header = FALSE, stringsAsFactors = FALSE)
rownames(autoimmune)<-autoimmune$V1

autoshared=matrix(0,dim(autoimmune)[1]*(dim(autoimmune)[1]-1)/2,1)
autoshared=data.table(autoshared)
write.table(autoshared, file = file.path(output_dir, "Autoimmune_sharedMHC.txt"), 
            quote=F,row.names=F,col.names = F,sep='\t')
autoshared <- read.table(file.path(output_dir, "Autoimmune_sharedMHC.txt"), 
                         header = FALSE, stringsAsFactors = FALSE)
autosharedindex=1
for (i in 1:(dim(autoimmune)[1]-1)){
  for (j in (i+1):dim(autoimmune)[1]){
    gene1 <- read.table(file.path(output_dir, paste0(autoimmune[i, 1], '_AssoGeneMHC.txt')), 
                        header=TRUE, stringsAsFactors = FALSE)
    gene2 <- read.table(file.path(output_dir, paste0(autoimmune[j, 1], '_AssoGeneMHC.txt')), 
                        header=TRUE, stringsAsFactors = FALSE)
    a=intersect(gene1$gene,gene2$gene)
    autoshared[autosharedindex,1]=length(a)
    rownames(autoshared)[autosharedindex]=paste(autoimmune[i,1],autoimmune[j,1],sep = '+')
    autosharedindex=autosharedindex+1
  }
}

write.table(autoshared, file = file.path(output_dir, "Autoimmune_shared_edgenumberMHC.txt"), 
            row.names = TRUE,col.names = FALSE,quote = FALSE)

for (i in 1:dim(autoimmune)[1]){
  dattmp=subset(datSigMHC,datSigMHC$StudyAccession==autoimmune[i,1])
  # for five autoimmune disorders, filtering the genes with low pred.perf.qvalue 
  # dattmp=subset(dattmp,dattmp$pred_perf_qval<=0.001)
  dattmp=subset(dattmp,dattmp$pva.qval<=0.5)
  dattmp=dattmp[,1:2]
  write.table(dattmp, 
              file = file.path(output_dir, paste0(autoimmune[i, 1], '_AssoGene_filterMHC.txt')), 
              quote = FALSE, row.names = FALSE, col.names = TRUE)
}


## filter genes with high effects (up or downregulated or both)

for (i in 1:(dim(autoimmune)[1])){
  gene1 <- read.table(file.path(output_dir, paste0(autoimmune[i, 1], '_AssoGene_filterMHC.txt')), 
                      header = TRUE, stringsAsFactors = FALSE)
  a=gene1$gene
  aa=paste('tmp1=',autoimmune[i,1],'[',autoimmune[i,1],'$gene %in% a,]',sep='')
  eval(parse(text = aa))
  # assogene=tmp1[abs(tmp1$zscore)>=mean(abs(tmp1$zscore)),]
  assogene =tmp1
  write.table(assogene[, 1:3], 
              file = file.path(output_dir, paste0(autoimmune[i, 1], '_SharedGeneHigheffectMHC.txt')), 
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}

autosharedindex=1
for (i in 1:(dim(autoimmune)[1]-1)){
  for (j in (i+1):dim(autoimmune)[1]){
    gene1 <- read.table(file.path(output_dir, paste0(autoimmune[i, 1], '_SharedGeneHigheffectMHC.txt')), 
                        header = TRUE, stringsAsFactors = FALSE)
    gene2 <- read.table(file.path(output_dir, paste0(autoimmune[j,1], '_SharedGeneHigheffectMHC.txt')), 
                        header = TRUE, stringsAsFactors = FALSE)
    a=intersect(gene1$gene,gene2$gene)
    write.table(a, 
                file = file.path(output_dir, paste0(autoimmune[i, 1],'+',autoimmune[j,1],'_SharedGeneHigheffectMHC.txt')), 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    autoshared[autosharedindex,1]=length(a)
    rownames(autoshared)[autosharedindex]=paste(autoimmune[i,1],autoimmune[j,1],sep = '+')
    autosharedindex=autosharedindex+1
  }
}


#a=unique(datSig1$gene)
#for (i in 1:(dim(autoimmune)[1]-1)){
#  for (j in (i+1):dim(autoimmune)[1]){
#    if (file.info(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffect.txt',sep=''))$size !=0){
#    gene1<-read.table(paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffect.txt',sep=''),header=FALSE,stringsAsFactors=FALSE)
#    a=intersect(gene1[,1],gene2$gene)
#    write.table(a,file=paste('/Users/wenzhang/Desktop/Trait_AGene/',autoimmune[i,1],'+',autoimmune[j,1],'_SharedGeneHigheffect.txt',sep=''),quote=FALSE,col.names = FALSE,row.names = FALSE)
#    autoshared[autosharedindex,1]=length(a)
#    rownames(autoshared)[autosharedindex]=paste(autoimmune[i,1],autoimmune[j,1],sep = '+')
#    autosharedindex=autosharedindex+1}

#  }
#}


autoimmuneGeneTrait=matrix(0,sum(autoshared)*2,1)
autoimmuneGeneTrait=data.table(autoimmuneGeneTrait)
write.table(autoimmuneGeneTrait, 
            file = file.path(output_dir, "Autoimmune_Gene_correlationMHC.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
autoimmuneGeneTrait <- read.table(file.path(output_dir,"Autoimmune_Gene_correlationMHC.txt"), 
                    header = FALSE, stringsAsFactors = FALSE)
autogenetraitindex=1

for (i in 1:(dim(autoimmune)[1] - 1)) {
  for (j in (i + 1):dim(autoimmune)[1]) {
    file_path <- file.path(output_dir, paste0(autoimmune[i, 1], '+', autoimmune[j, 1], '_SharedGeneHigheffectMHC.txt'))
    
    # 检查文件是否存在且大小不为零
    if (file.exists(file_path) && file.info(file_path)$size != 0) {
      gene1 <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
      
      # 检查 gene1 是否有数据行
      if (dim(gene1)[1] > 0) {
        autoimmuneGeneTrait[autogenetraitindex:(dim(gene1)[1] + autogenetraitindex - 1), 1] <- autoimmune[i, 1]
        autoimmuneGeneTrait[autogenetraitindex:(dim(gene1)[1] + autogenetraitindex - 1), 2] <- gene1$V1
        autoimmuneGeneTrait[(dim(gene1)[1] + autogenetraitindex):(2 * dim(gene1)[1] + autogenetraitindex - 1), 1] <- autoimmune[j, 1]
        autoimmuneGeneTrait[(dim(gene1)[1] + autogenetraitindex):(2 * dim(gene1)[1] + autogenetraitindex - 1), 2] <- gene1$V1
        autogenetraitindex <- autogenetraitindex + 2 * dim(gene1)[1]
      }
    }
  }
}

# set node size
autoimmuneGeneTrait[,3]=autoimmune[autoimmuneGeneTrait$V1,2]
# set trait color (corresponds to yellow)
autoimmuneGeneTrait[,4]=4
# set gene color (corresponds to grey)
autoimmuneGeneTrait[,5]=5
# set gene size (should be the same)
autoimmuneGeneTrait[,6]=2.5

#############################################################################
allnumber=0
for(i in 1:dim(autoimmuneGeneTrait)[1]){
  aa=(datSigMHC[datSigMHC$StudyAccession==autoimmuneGeneTrait[i,1],])
  aa=aa[aa$gene==autoimmuneGeneTrait[i,2],]
  # set edgewidth: number of traits (tissues via which) the gene is highly correlated with
  autoimmuneGeneTrait[i,7]=dim(aa)[1]
  # raw evaluation of up/down regulations
  upordow=(aa$zscore>=0)
  allnumber=0
  for (j in 1:length(upordow))
  {
    if (isTRUE(upordow[j])) allnumber=allnumber+1
  }
  # all positively correlated
  if (allnumber==length(upordow)) autoimmuneGeneTrait[i,8]=1
  # there are tissues where the gene is negatively correlated with the trait 
  if (allnumber<length(upordow)) {
    autoimmuneGeneTrait[i,8]=3
    # positively correlated cases are less than half --> negatively correlated
    if (allnumber<length(upordow)/2) autoimmuneGeneTrait[i,8]=2
  }
  
}

for(i in 1:dim(autoimmuneGeneTrait)[1]){
  aa=datSigMHC[datSigMHC$gene==autoimmuneGeneTrait[i,2],]
  autoimmuneGeneTrait[i,2]=aa[1,2]
}

# advanced evaluations of up/down regulations
autoimmuneGeneReg=autoimmuneGeneTrait[,1:2]
autoimmuneGeneReg[,3:52]=0
colnames(autoimmuneGeneReg)<-c('StudyAccession','gene',unique(datSigMHC1$tissue))
for(i in 1:dim(autoimmuneGeneReg)[1]){
  aa=(datSigMHC[datSigMHC$StudyAccession==autoimmuneGeneReg[i,1],])
  aa=aa[aa$gene_name==autoimmuneGeneReg[i,2],]
  for (j in 3:16){
    tmp=aa[aa$tissue==colnames(autoimmuneGeneReg)[j],]
    if (dim(tmp)[1]>0){
      if (tmp$zscore<0){
        autoimmuneGeneReg[i,j]=1
      }
      if (tmp$zscore>0){
        autoimmuneGeneReg[i,j]=2
      }
    }
  }
  
  # the number counts of no, down and up regulations:
  a0=0
  a1=0
  a2=0
  for (j in 3:52){
    if (autoimmuneGeneReg[i,j]==0) a0=a0+1
    if (autoimmuneGeneReg[i,j]==1) a1=a1+1
    if (autoimmuneGeneReg[i,j]==2) a2=a2+1
  }
  if (a0==49) autoimmuneGeneReg[i,52]=0
  if (a0!=49){
    if (a1==a2) autoimmuneGeneReg[i,52]=0
    if (a1>a2) autoimmuneGeneReg[i,52]=1
    if (a1<a2)  autoimmuneGeneReg[i,52]=2
  }
  
  
  # # set edgewidth: number of traits (tissues via which) the gene is highly correlated with
  #  autoimmuneGeneTrait[i,7]=dim(aa)[1]
  #  # raw evaluation of up/down regulations
  # upordow=(aa$zscore>=0)
  #  allnumber=0
  #  for (j in 1:length(upordow))
  #  {
  #   if (isTRUE(upordow[j])) allnumber=allnumber+1
  #  }
  # # all positively correlated
  #  if (allnumber==length(upordow)) autoimmuneGeneTrait[i,8]=1
  #  # there are tissues where the gene is negatively correlated with the trait 
  #  if (allnumber<length(upordow)) {
  #    autoimmuneGeneTrait[i,8]=3
  #    # positively correlated cases are less than half --> negatively correlated
  #    if (allnumber<length(upordow)/2) autoimmuneGeneTrait[i,8]=2
  #  }
  
}

colnames(autoimmuneGeneReg)[52]<-'Regulation(1down/2up)'

autoimmuneGeneTrait[,8]<-autoimmuneGeneReg$`Regulation(1down/2up)`

autoimmuneGeneTrait=autoimmuneGeneTrait[,-5]

colnames(autoimmuneGeneTrait)<-c('StudyAccession','gene','node1Size','Node1color','Node2Color/Size','edgeWidth','Regulation(1down/2up)')

autoimmuneGeneTrait<-autoimmuneGeneTrait[autoimmuneGeneTrait$`Regulation(1down/2up)`!=0,]

write.table(autoimmuneGeneTrait, 
            file = file.path(output_dir, "Autoimmune_Trait_geneMHC.txt"), 
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')