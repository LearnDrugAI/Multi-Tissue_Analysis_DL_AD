#### for each trait, get number of genes that are up/down regulated
# 加载必要的库
library(qvalue)
library(RColorBrewer)
library(WGCNA)

# 读取数据
All_data_clean <- read.csv('All_data_clean.csv', header = TRUE, stringsAsFactors = FALSE, sep = ",")
datMHCmi <- All_data_clean

# 计算 q-values 并筛选显著结果
qobjwo <- qvalue(datMHCmi$pvalue, fdr.level = 0.01)
datMHCmi$pva.qval <- qobjwo$qvalues
datMHCmi$isSignificant <- qobjwo$significant
datMHCmi <- datMHCmi[datMHCmi$pva.qval <= 0.01, ]

# 选择显著的研究项目
selected_studies <- c("GCST90041866", "GCST90041879", "GCST90042665", "GCST90042678", "GCST90042691")

# 根据逻辑条件筛选出符合条件的行
datMHCmi <- datMHCmi[datMHCmi$StudyAccession %in% selected_studies, ]

# 生成显著性计数表
datTabMHC <- table(datMHCmi[datMHCmi$isSignificant, c("tissue", "StudyAccession")])

datSigMHC <- datMHCmi[datMHCmi$isSignificant, ]

traitupdownGene=matrix(0,5,3)

rownames(traitupdownGene)<- c('GCST90041866','GCST90041879','GCST90042665','GCST90042678','GCST90042691')
# second column: number of up regulated genes
# third column: number of down regulated genes
# fourth column: number of ambiguous regulated genes

for ( i in 1:dim(traitupdownGene)[1]){
  dattmp=datSigMHC[datSigMHC$StudyAccession==rownames(traitupdownGene)[i],]
  # genes that associated with the trait
  onlyassogene=unique(dattmp$gene)
  for(j in 1:length(onlyassogene)){
    # for each gene: judge if it is up/down/ambiguous expressional changed during the trait
    # get all tissues that the gene is associated with the trait
    setgene=dattmp[dattmp$gene==onlyassogene[j],]
    
    unordow=matrix(0,1,49)
    colnames(unordow)<-unique(datSigMHC$tissue)
    upchange=0
    downchange=0
    for (k in 1:49) {
      tmp=setgene[setgene$tissue==colnames(unordow)[k],]
      print(setgene[setgene$tissue==colnames(unordow)[49],])
      if (dim(tmp)[1]>0){
        if (tmp$zscore<0){ # down: set value 1
          unordow[k]=1
        }
        if (tmp$zscore>0){ # up: set value 2
          unordow[k]=2
        }
      }
      
      if (unordow[k]!=0) {
        if (unordow[k]==2) upchange=upchange+1
        if (unordow[k]==1) downchange=downchange+1    }
      
    }
    if (upchange>downchange) traitupdownGene[i,1]=traitupdownGene[i,1]+1
    if(upchange<downchange)  traitupdownGene[i,2]=traitupdownGene[i,2]+1
    if(upchange==downchange)  traitupdownGene[i,3]=traitupdownGene[i,3]+1
    
  }
}


# indcategory=rep(0,5)
# for (i in 1:length(Traitname[,1])){
#   for (j in 1:length(rownames(traitupdownGene))){
#     if (rownames(traitupdownGene)[j]==Traitname[i,1])
#       indcategory[i]=j
#   }
# }

# traitupdownGene=traitupdownGene[indcategory,]
write.table(traitupdownGene,'Result/TraitGeneNumber.txt',quote=FALSE,col.names = FALSE)
traitupdownGene=read.table('Result/TraitGeneNumber.txt',header=FALSE,stringsAsFactors = FALSE)
rownames(traitupdownGene)<-traitupdownGene$V1
traitupdownGene[,5]=Traitname[rownames(traitupdownGene),2]
traitupdownGene=traitupdownGene[,-1]

# Trait<-traitupdownGene[,4]
# Trait<-c(rep(Trait,each=3))
# Regulation<-c(rep(c("Up", "Down", "Ambiguous"), times = 58))
# Frequency<-rep(0,58*3)
# for(i in 1:58){
#   for (j in 1:3 ){
#     Frequency[(i-1)*3+j]=traitupdownGene[i,j]
#   }
# }
# Data <- data.frame(Trait, Regulation, Frequency)
# 
# mpdf('barchar')
#  ggplot(Data, aes(x = Trait, y = Frequency, fill = Regulation, label = Frequency)) +
#  geom_bar(stat = "identity") +
#   geom_text(size = 3, position = position_stack(vjust = 90))
#  dev.off()

install.packages("plotly")
library(plotly)

#  data <- data.frame(traitupdownGene[,4], traitupdownGene[,1], traitupdownGene[,2],traitupdownGene[,3])
#  colnames(data)<-c('Trait','Up','Down','Ambiguous')
#  #The default order will be alphabetized unless specified as below:
#  data$Trait <- factor(data$Trait, levels = data[["Trait"]])
#  mpdf("BarPlot", width=2+58*0.30, height=4+length(60)*0.30)
#  mar.default = c(6,5,5,3) + 0.1
#  par(mar = mar.default + c(10, 12, 0, 0)) #c(bottom, left, top, right)
# plot_ly(data, x = ~Trait, y = ~Down, type = 'bar', name = 'Downregulation', marker = list(color = 'rgb(49,130,189)')) %>%
#   add_trace(y = ~Ambiguous, name = 'Ambiguous', marker = list(color = 'rgb(4,154,4)')) %>%
#   add_trace(y = ~Up,  name = 'Upregulation', marker = list(color = 'rgb(204,0,14)')) %>%
#      layout(xaxis = list(title = "", tickangle = -65),
#           yaxis = list(title = "Number of correlated genes"),font=t,
#           margin = list(b = 225),
#           barmode = 'group')
#  dev.off()
# 创建一个数据框
data <- data.frame(Trait = rownames(traitupdownGene), 
                   Up = traitupdownGene[,1], 
                   Down = traitupdownGene[,2],
                   Ambiguous = traitupdownGene[,3])

# 设置数据框列名
colnames(data) <- c('Trait','Up','Down','Ambiguous')

# 设置Trait列为因子，以保持顺序
data$Trait <- factor(data$Trait, levels = data[["Trait"]])

# 创建交互式条形图
plot <- plot_ly(data, x = ~Trait, y = ~Down, type = 'bar', name = 'Downregulation', marker = list(color = 'rgb(49,130,189)')) %>%
  add_trace(y = ~Ambiguous, name = 'Ambiguous', marker = list(color = 'rgb(4,154,4)')) %>%
  add_trace(y = ~Up,  name = 'Upregulation', marker = list(color = 'rgb(204,0,14)')) %>%
  layout(xaxis = list(title = "", tickangle = -65),
         yaxis = list(title = "Number of correlated genes",
                      titlefont = list(family = "Times New Roman", size = 14)), # 设置y轴标签的字体和大小
         margin = list(t = 100), # 调整上边距和下边距，将图移到中间
         barmode = 'group',
         font = list(family = "Times New Roman")) # 设置图中所有字体的样式为 Times New Roman

# 使用pdf函数创建PDF文件
pdf("BarPlot2.pdf", width = 8.5, height = 11)  # 设置页面宽度和高度

# 打印交互式条形图
print(plot)

# 关闭PDF设备
dev.off()

# 对于每个特征
for (i in 1:dim(traitupdownGene)[1]) {
  dattmp <- datSigMHC[datSigMHC$StudyAccession == rownames(traitupdownGene)[i],]
  
  # 提取与特征相关的基因和基因名称
  onlyassogene <- unique(dattmp$gene)
  gene_names <- unique(dattmp$gene_name)
  
  # 初始化上调和下调基因向量
  up_genes_for_trait <- data.frame(gene = character(0), gene_name = character(0))
  down_genes_for_trait <- data.frame(gene = character(0), gene_name = character(0))
  
  # 对于每个基因
  for (j in 1:length(onlyassogene)) {
    setgene <- dattmp[dattmp$gene == onlyassogene[j],]
    
    # 初始化上调和下调标志
    upchange <- 0
    downchange <- 0
    
    # 对于每个组织
    for (k in 1:49) {
      tmp <- setgene[setgene$tissue == unique(datSigMHC$tissue)[k],]
      
      # 如果存在数据
      if (dim(tmp)[1] > 0) {
        if (tmp$zscore < 0) { # 下调
          downchange <- 1
        }
        if (tmp$zscore > 0) { # 上调
          upchange <- 1
        }
      }
    }
    
    # 如果存在上调或下调，则将基因及其名称添加到相应的列表中
    if (upchange == 1) {
      up_genes_for_trait <- rbind(up_genes_for_trait, data.frame(gene = onlyassogene[j], gene_name = gene_names[j]))
    }
    if (downchange == 1) {
      down_genes_for_trait <- rbind(down_genes_for_trait, data.frame(gene = onlyassogene[j], gene_name = gene_names[j]))
    }
  }
  
  # 将上调和下调基因及其名称列表保存到全局列表中
  up_genes[[i]] <- up_genes_for_trait
  down_genes[[i]] <- down_genes_for_trait
}

# 将上调和下调基因及其名称列表保存为文件
for (i in 1:length(up_genes)) {
  write.table(up_genes[[i]], paste0("Result/updowngene/", rownames(traitupdownGene)[i], "_UpGenes.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  write.table(down_genes[[i]], paste0("Result/updowngene/", rownames(traitupdownGene)[i], "_DownGenes.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}