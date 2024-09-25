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

TRAIT=unique(datMHC$StudyAccession)
TISSUE=unique(datMHC$tissue)
GENE=unique(datMHC$gene)

# indcategory=rep(0,61)
# for (i in 1:length(TRAIT)) {
#   for (j in 1:length(colnames(datTabMHC))) {
#     if (colnames(datTabMHC)[j] == TRAIT[i]) {
#       indcategory[i] = j
#       break  # 找到匹配项后，退出内循环
#     }
#   }
# }



Tissuename <- data.frame(TISSUE)
rownames(Tissuename)<-Tissuename$'TISSUE'

Traitname <- data.frame(TRAIT)
rownames(Traitname)<-Traitname$'TRAIT'
# all of the numbers
datSigMHC=datMHC[datMHC$isSignificant,]
#datSigMHC1=datSigMHC[abs(datSigMHC$zscore)>min(abs(datSigMHC$zscore)),]

#allTabMHC1=table(datSigMHC1[,c("tissue","trait")])
sigTabMHC=table(datSigMHC[,c("tissue","StudyAccession")])
allTabMHC=table(datMHC[,c("tissue","StudyAccession")])#trait改成了StudyAccession

# normalized count numbers
# 使用intersect函数获取两个表中共有的行和列
common_rows <- intersect(rownames(sigTabMHC), rownames(allTabMHC))
common_cols <- intersect(colnames(sigTabMHC), colnames(allTabMHC))

# 使用公共的行和列进行子集选择
sigTabMHC_subset <- sigTabMHC[common_rows, common_cols]
allTabMHC_subset <- allTabMHC[common_rows, common_cols]

# 计算标准化的显著计数
NORMALIZED_SIGNIFICANT_COUNTS_MHC <- sigTabMHC_subset / allTabMHC_subset


# scale the normalized count
matrix_normMHC=NORMALIZED_SIGNIFICANT_COUNTS_MHC
for(i in 1:49){
  for (j in 1:18){
    matrix_normMHC[i,j]=(NORMALIZED_SIGNIFICANT_COUNTS_MHC[i,j]-mean(NORMALIZED_SIGNIFICANT_COUNTS_MHC[,j]))/sd(NORMALIZED_SIGNIFICANT_COUNTS_MHC[,j])
  }
}#将14修改为49，将58修改为18


#trait改成了StudyAccession



traitGenecountMHC=matrix(0,length(TRAIT),1)
tissueGenecountMHC=matrix(0,length(TISSUE),1)
geneTraitcountMHC=matrix(0,length(GENE),1)


rownames(traitGenecountMHC)<-TRAIT
rownames(tissueGenecountMHC)<-TISSUE
rownames(geneTraitcountMHC)<-GENE


for (i in 1:length(GENE)){
  dattmp=subset(datSigMHC,datSigMHC$gene==GENE[i])
  geneTraitcountMHC[i]=length(unique(dattmp$StudyAccession))
}

for (i in 1:length(TRAIT)){
  dattmp=subset(datSigMHC,datSigMHC$StudyAccession==TRAIT[i])
  traitGenecountMHC[i]=length(unique(dattmp$gene))
}

for (i in 1:length(TISSUE)){
  dattmp=subset(datSigMHC,datSigMHC$tissue==TISSUE[i])
  tissueGenecountMHC[i]=length(unique(dattmp$gene))
}



geneTraitcoun=geneTraitcountMHC
colnames(geneTraitcoun)<-'Traitnumber'
geneTraitcoun=geneTraitcoun[order(geneTraitcoun,decreasing = TRUE)]



colNamMHC<-colnames(matrix_normMHC)
colNam1MHC=Traitname[colNamMHC,1]
colNam1MHC<-paste(colNam1MHC," (",traitGenecountMHC[colNamMHC,],") ",sep='')
rowNamMHC<-rownames(matrix_normMHC)
rowNam1MHC=Tissuename[rowNamMHC,1]
rowNam1MHC<-paste(rowNam1MHC," (",tissueGenecountMHC[rowNamMHC,],") ",sep='')

# matrix_normMHC=matrix_normMHC[,indcategory]
# datTabMHC=datTabMHC[,indcategory]


pdf("SIGNIFICANT_COUNTS.NORMALIZED_RATIO_countNumber1-MHC_test.pdf", width=3+length(unique(datMHC$StudyAccession))*0.30, height=4.5+length(unique(datMHC$tissue))*0.30)
mar.default = c(6,5,5,3) + 0.1
par(mar = mar.default + c(11, 13, 0, 0)) #c(bottom, left, top, right)
labeledHeatmap(Matrix = apply(matrix_normMHC, 2, scale),
               xLabels = colNam1MHC,
               yLabels = rowNam1MHC,
               colorLabels = F,
               colors = myPalette2way(1000),
               textMatrix = datTabMHC,
               setStdMargins = FALSE,
               cex.text = 0.65, 
               zlim = c(-3.5,3.5),
               naColor = "white",
               main = "")
dev.off()


