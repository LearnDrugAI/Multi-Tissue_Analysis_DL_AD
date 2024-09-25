library(biomaRt)
library(data.table)

# 加载基因表达数据
mtx.rna <-read.table(file = "geneExpression.csv", header = TRUE, sep = ",",row.names = 1, check.names = FALSE)
# 如果行名是基因名称，则使用 rownames 提取基因名称
genes <- rownames(mtx.rna)

# 查看基因名称
head(genes)

# 连接到 Ensembl 数据库
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# 获取基因注释信息
gene_annot <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "external_gene_name",
  values = genes,
  mart = mart
)

# 重命名列以符合要求
colnames(gene_annot) <- c("chr", "start", "end", "gene_id", "gene_name", "gene_type")

# 选择并重新排列列的顺序
gene_annot <- gene_annot[, c("chr", "gene_id", "gene_name", "start", "end", "gene_type")]

# 查看结果
head(gene_annot)

# 保存为文件
write.table(gene_annot, "gene_annot.parsed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

######################filt##############################################
# 加载数据
gene_annot <- fread("gene_annot.parsed.txt", header = TRUE, sep = "\t")

# 查看染色体名称
unique(gene_annot$chr)

# 过滤掉非常染色体的数据
filtered_gene_annot <- gene_annot[gene_annot$chr %in% as.character(1:22), ]

# 保存过滤后的数据
write.table(filtered_gene_annot, "gene_annot.parsed.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#test = read.table(file = "genotype.chr22.txt", header = TRUE, sep = "\t",nrows = 2)

