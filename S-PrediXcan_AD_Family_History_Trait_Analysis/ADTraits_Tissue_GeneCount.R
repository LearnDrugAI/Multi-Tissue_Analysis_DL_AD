library(qvalue)
library(RColorBrewer)
library(WGCNA)
library(data.table)

outDir='/Result'

All_data_clean <- fread('All_data_clean.csv', header = TRUE, stringsAsFactors = FALSE, sep = ",")
# data that reduce SNPs in MHC region
datMHC=All_data_clean


myPalette = colorRampPalette(brewer.pal(9, "Greens"), space="Lab")
myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")
my_palette=colorRampPalette(c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B"),space="Lab")


## for removing MHC results
#match=is.na(datMHC$pvalue)
#datMHC=datMHC[!match,]
qobjwo <- qvalue(datMHC$pvalue, fdr.level = 0.01)
datMHC$pva.qval=qobjwo$qvalues
datMHC$isSignificant=qobjwo$significant
datMHC=datMHC[datMHC$pva.qval<=0.01,] #pred_perf_qval改成了pva.qval


datTabMHC=table(datMHC[datMHC$isSignificant,c("tissue","StudyAccession")]) #trait改成了StudyAccession

# 创建一个包含目标 StudyAccession 值的向量
selected_study_accessions <- c('GCST90041866', 'GCST90041879', 'GCST90042665', 'GCST90042678', 'GCST90042691')

# 根据逻辑条件筛选出符合条件的行
selected_rows <- datSigMHC[datSigMHC$StudyAccession %in% selected_study_accessions, ]
#write.csv(selected_rows, "datBra_5important.csv", row.names = FALSE)


TRAIT=unique(selected_rows$StudyAccession)
TISSUE=unique(selected_rows$tissue)
GENE=unique(selected_rows$gene)
#trait改成了StudyAccession

# traitGenecount=matrix(0,length(TRAIT),1)
# tissueGenecount=matrix(0,length(TISSUE),1)
# geneTraitcount=matrix(0,length(GENE),1)

traitGenecountMHC=matrix(0,length(TRAIT),1)
tissueGenecountMHC=matrix(0,length(TISSUE),1)
geneTraitcountMHC=matrix(0,length(GENE),1)


# rownames(traitGenecount)<-TRAIT
# rownames(tissueGenecount)<-TISSUE
# rownames(geneTraitcount)<-GENE

rownames(traitGenecountMHC)<-TRAIT
rownames(tissueGenecountMHC)<-TISSUE
rownames(geneTraitcountMHC)<-GENE

for (i in 1:length(GENE)){
  dattmp=subset(selected_rows,selected_rows$gene==GENE[i])
  geneTraitcountMHC[i]=length(unique(dattmp$StudyAccession))
}

for (i in 1:length(TRAIT)){
  dattmp=subset(selected_rows,selected_rows$StudyAccession==TRAIT[i])
  traitGenecountMHC[i]=length(unique(dattmp$gene))
}

for (i in 1:length(TISSUE)){
  dattmp=subset(selected_rows,selected_rows$tissue==TISSUE[i])
  tissueGenecountMHC[i]=length(unique(dattmp$gene))
}

# 导入 ggplot2 库用于绘图
library(ggplot2)

# 将矩阵转换为数据框
trait_genecount_df <- data.frame(StudyAccession = TRAIT, GeneCount = traitGenecountMHC)
tissue_genecount_df <- data.frame(Tissue = TISSUE, GeneCount = tissueGenecountMHC)


# 画图：Tissue vs GeneCount
pdf("Tissue_GeneCount_Distribution.pdf", width = 10, height = 6)  # 调整 PDF 尺寸
ggplot(tissue_genecount_df, aes(x = reorder(Tissue, GeneCount), y = GeneCount)) +
  geom_bar(stat = "identity", fill = "lightgreen", width = 0.5) +
  labs(title = "Number of Genes per Tissue",
       x = "Tissue",
       y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))  # 缩小字体并旋转 x 轴标签
dev.off()


# 统一参数
common_theme <- theme_classic() +  # 使用theme_classic去掉背景的灰色
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # 标题字体放大
    axis.title.x = element_text(size = 14),  # x轴标题放大
    axis.title.y = element_text(size = 14),  # y轴标题放大
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # x轴标签旋转并放大
    axis.text.y = element_text(size = 12)  # y轴标签放大
  )

# 画图：StudyAccession vs GeneCount
# 绘制StudyAccession的图表
pdf("Trait_GeneCount_Distribution.pdf", width = 10, height = 6)  # 调整PDF尺寸
ggplot(trait_genecount_df, aes(x = StudyAccession, y = GeneCount)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  geom_text(aes(label = GeneCount), vjust = -0.5, size = 4) +  # 添加具体数量标签并放大
  labs(title = "Number of Genes per StudyAccession",
       x = "StudyAccession",
       y = "Number of Genes") +
  common_theme  # 统一的主题
dev.off()


